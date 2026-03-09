/*
 * metal_m2l.mm
 * Objective-C++ implementation of Metal-accelerated M2L transformation
 *
 * Precision: Float + Kahan summation (~7 digits with error compensation)
 * For full double precision (~15 digits), use CPU-based M2L instead.
 *
 * Compile with Apple Clang:
 *   clang++ -c metal_m2l.mm -o metal_m2l.o -std=c++17 -fobjc-arc
 *   clang++ -shared -o libmetal_m2l.dylib metal_m2l.o -framework Metal -framework Foundation
 */

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#include "metal_m2l.h"
#include <cstring>
#include <vector>

/* Context structure */
struct MetalM2LContext {
    id<MTLDevice> device;
    id<MTLCommandQueue> commandQueue;
    id<MTLComputePipelineState> pipelineFloatKahan;
    id<MTLComputePipelineState> pipelineFloatKahanTG;
    id<MTLLibrary> library;
    int threadgroupMode;  /* 0=off, 1=TG */

    /* Static data buffers (set up once when FMM structure is initialized) */
    id<MTLBuffer> rowOffsetBuffer;
    id<MTLBuffer> rowLenBuffer;

    /* Float precision coefficient buffers */
    id<MTLBuffer> coeffReBuffer;
    id<MTLBuffer> coeffImBuffer;

    /* Source index buffer */
    id<MTLBuffer> termSrcIdxBuffer;

    /* Source MM buffers - float precision */
    id<MTLBuffer> srcMmRe0Buffer;
    id<MTLBuffer> srcMmIm0Buffer;
    id<MTLBuffer> srcMmRe1Buffer;
    id<MTLBuffer> srcMmIm1Buffer;

    /* Result buffer */
    id<MTLBuffer> resultBuffer;

    /* For async execution */
    id<MTLCommandBuffer> pendingCommandBuffer;

    /* Dimensions */
    int32_t numBuckets;
    int32_t numCoeffs;
    int32_t numRows;
    int32_t numTerms;
    int32_t totalMmSize;  /* numBuckets * numCoeffs */

    bool setupDone;

    char deviceName[256];
};

extern "C" {

int metal_m2l_is_available(void) {
    @autoreleasepool {
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        return (device != nil) ? 1 : 0;
    }
}

MetalM2LContext* metal_m2l_init(int threadgroup_mode) {
    @autoreleasepool {
        MetalM2LContext* ctx = new MetalM2LContext();
        memset(ctx, 0, sizeof(MetalM2LContext));

        /* Get the default Metal device */
        ctx->device = MTLCreateSystemDefaultDevice();
        if (ctx->device == nil) {
            delete ctx;
            return nullptr;
        }

        /* Store device name */
        const char* name = [[ctx->device name] UTF8String];
        strncpy(ctx->deviceName, name, sizeof(ctx->deviceName) - 1);
        ctx->deviceName[sizeof(ctx->deviceName) - 1] = '\0';

        /* Create command queue */
        ctx->commandQueue = [ctx->device newCommandQueue];
        if (ctx->commandQueue == nil) {
            delete ctx;
            return nullptr;
        }

        /* Load the Metal shader library */
        NSError* error = nil;

        /* Try to load from metallib first */
        NSString* libPath = @"metal_m2l.metallib";
        NSURL* libURL = [NSURL fileURLWithPath:libPath];
        ctx->library = [ctx->device newLibraryWithURL:libURL error:&error];

        if (ctx->library == nil) {
            /* Fallback: compile from source at runtime */
            NSFileManager* fm = [NSFileManager defaultManager];
            NSArray* searchPaths = @[
                @".",
                @"./metal_m2l",
                [[NSBundle mainBundle] bundlePath],
                @"/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_bem/metal_m2l"
            ];

            NSString* foundPath = nil;
            for (NSString* dir in searchPaths) {
                NSString* fullPath = [dir stringByAppendingPathComponent:@"metal_m2l.metal"];
                if ([fm fileExistsAtPath:fullPath]) {
                    foundPath = fullPath;
                    break;
                }
            }

            if (foundPath != nil) {
                NSString* source = [NSString stringWithContentsOfFile:foundPath
                                                             encoding:NSUTF8StringEncoding
                                                                error:&error];
                if (source != nil) {
                    MTLCompileOptions* options = [[MTLCompileOptions alloc] init];
                    options.mathMode = MTLMathModeSafe;  // Preserve Kahan summation precision
                    ctx->library = [ctx->device newLibraryWithSource:source
                                                             options:options
                                                               error:&error];
                }
            }
        }

        if (ctx->library == nil) {
            NSLog(@"[Metal M2L] Failed to load Metal library: %@", error);
            delete ctx;
            return nullptr;
        }

        /* Create float+Kahan pipeline */
        id<MTLFunction> floatKahanFunc = [ctx->library newFunctionWithName:@"m2l_float_kahan"];
        if (floatKahanFunc != nil) {
            ctx->pipelineFloatKahan = [ctx->device newComputePipelineStateWithFunction:floatKahanFunc error:&error];
            if (ctx->pipelineFloatKahan == nil) {
                NSLog(@"[Metal M2L] Warning: Failed to create float+Kahan pipeline: %@", error);
            }
        }

        /* Create float+Kahan threadgroup pipeline */
        id<MTLFunction> floatKahanTGFunc = [ctx->library newFunctionWithName:@"m2l_float_kahan_tg"];
        if (floatKahanTGFunc != nil) {
            ctx->pipelineFloatKahanTG = [ctx->device newComputePipelineStateWithFunction:floatKahanTGFunc error:&error];
            if (ctx->pipelineFloatKahanTG == nil) {
                NSLog(@"[Metal M2L] Warning: Failed to create float+Kahan TG pipeline: %@", error);
            }
        }

        if (ctx->pipelineFloatKahan == nil && ctx->pipelineFloatKahanTG == nil) {
            NSLog(@"[Metal M2L] Error: No valid pipeline created");
            delete ctx;
            return nullptr;
        }

        ctx->threadgroupMode = 0;
        ctx->setupDone = false;

        /* Set threadgroup mode from parameter (read from settings.json) */
        if (threadgroup_mode >= 1 && ctx->pipelineFloatKahanTG != nil) {
            ctx->threadgroupMode = 1;
            NSLog(@"[Metal M2L] Threadgroup kernel enabled (TG=1)");
        } else if (threadgroup_mode != 0) {
            NSLog(@"[Metal M2L] Warning: TG kernel mode %d requested but unavailable", threadgroup_mode);
        }

        NSLog(@"[Metal M2L] Initialized on: %s (TG=%d)", ctx->deviceName, ctx->threadgroupMode);
        return ctx;
    }
}

int metal_m2l_setup(
    MetalM2LContext* ctx,
    int32_t num_buckets,
    int32_t num_coeffs,
    int32_t num_rows,
    int32_t num_terms,
    const int32_t* row_offset,
    const int32_t* row_len,
    const float* term_coeff_re,
    const float* term_coeff_im,
    const int32_t* term_src_idx)
{
    if (ctx == nullptr || (ctx->pipelineFloatKahan == nil && ctx->pipelineFloatKahanTG == nil)) {
        return -1;
    }

    @autoreleasepool {
        ctx->numBuckets = num_buckets;
        ctx->numCoeffs = num_coeffs;
        ctx->numRows = num_rows;
        ctx->numTerms = num_terms;
        ctx->totalMmSize = num_buckets * num_coeffs;

        /* Create row offset/len buffers */
        ctx->rowOffsetBuffer = [ctx->device newBufferWithBytes:row_offset
                                                        length:num_rows * sizeof(int32_t)
                                                       options:MTLResourceStorageModeShared];
        ctx->rowLenBuffer = [ctx->device newBufferWithBytes:row_len
                                                     length:num_rows * sizeof(int32_t)
                                                    options:MTLResourceStorageModeShared];

        /* Create coefficient buffers */
        ctx->coeffReBuffer = [ctx->device newBufferWithBytes:term_coeff_re
                                                      length:num_terms * sizeof(float)
                                                     options:MTLResourceStorageModeShared];
        ctx->coeffImBuffer = [ctx->device newBufferWithBytes:term_coeff_im
                                                      length:num_terms * sizeof(float)
                                                     options:MTLResourceStorageModeShared];

        /* Create source index buffer */
        ctx->termSrcIdxBuffer = [ctx->device newBufferWithBytes:term_src_idx
                                                         length:num_terms * sizeof(int32_t)
                                                        options:MTLResourceStorageModeShared];

        /* Create source MM buffers (will be filled each iteration) */
        size_t mmSize = ctx->totalMmSize * sizeof(float);
        ctx->srcMmRe0Buffer = [ctx->device newBufferWithLength:mmSize options:MTLResourceStorageModeShared];
        ctx->srcMmIm0Buffer = [ctx->device newBufferWithLength:mmSize options:MTLResourceStorageModeShared];
        ctx->srcMmRe1Buffer = [ctx->device newBufferWithLength:mmSize options:MTLResourceStorageModeShared];
        ctx->srcMmIm1Buffer = [ctx->device newBufferWithLength:mmSize options:MTLResourceStorageModeShared];

        /* Create result buffer */
        ctx->resultBuffer = [ctx->device newBufferWithLength:num_rows * 4 * sizeof(float)
                                                     options:MTLResourceStorageModeShared];

        ctx->setupDone = true;
        NSLog(@"[Metal M2L] Setup: %d buckets, %d coeffs, %d rows, %d terms [float+Kahan]",
              num_buckets, num_coeffs, num_rows, num_terms);
        return 0;
    }
}

/* Buffer accessors */
float* metal_m2l_get_src_mm_re0_buffer(MetalM2LContext* ctx) {
    return ctx && ctx->srcMmRe0Buffer ? (float*)[ctx->srcMmRe0Buffer contents] : nullptr;
}

float* metal_m2l_get_src_mm_im0_buffer(MetalM2LContext* ctx) {
    return ctx && ctx->srcMmIm0Buffer ? (float*)[ctx->srcMmIm0Buffer contents] : nullptr;
}

float* metal_m2l_get_src_mm_re1_buffer(MetalM2LContext* ctx) {
    return ctx && ctx->srcMmRe1Buffer ? (float*)[ctx->srcMmRe1Buffer contents] : nullptr;
}

float* metal_m2l_get_src_mm_im1_buffer(MetalM2LContext* ctx) {
    return ctx && ctx->srcMmIm1Buffer ? (float*)[ctx->srcMmIm1Buffer contents] : nullptr;
}

float* metal_m2l_get_result_buffer(MetalM2LContext* ctx) {
    return ctx && ctx->resultBuffer ? (float*)[ctx->resultBuffer contents] : nullptr;
}

int metal_m2l_compute_no_copy(MetalM2LContext* ctx) {
    if (ctx == nullptr || !ctx->setupDone) return -1;
    if (ctx->numRows <= 0) return 0;

    @autoreleasepool {
        id<MTLCommandBuffer> commandBuffer = [ctx->commandQueue commandBuffer];
        id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];
        id<MTLComputePipelineState> pipeline = nil;
        bool useThreadgroupKernel = false;

        if (ctx->threadgroupMode == 1 && ctx->pipelineFloatKahanTG != nil) {
            /* Float precision with threadgroup reduction */
            pipeline = ctx->pipelineFloatKahanTG;
            useThreadgroupKernel = true;
        } else if (ctx->pipelineFloatKahan != nil) {
            /* Float precision with Kahan summation */
            pipeline = ctx->pipelineFloatKahan;
        } else {
            [encoder endEncoding];
            return -1;
        }

        [encoder setComputePipelineState:pipeline];
        [encoder setBuffer:ctx->srcMmRe0Buffer offset:0 atIndex:0];
        [encoder setBuffer:ctx->srcMmIm0Buffer offset:0 atIndex:1];
        [encoder setBuffer:ctx->srcMmRe1Buffer offset:0 atIndex:2];
        [encoder setBuffer:ctx->srcMmIm1Buffer offset:0 atIndex:3];
        [encoder setBuffer:ctx->coeffReBuffer offset:0 atIndex:4];
        [encoder setBuffer:ctx->coeffImBuffer offset:0 atIndex:5];
        [encoder setBuffer:ctx->termSrcIdxBuffer offset:0 atIndex:6];
        [encoder setBuffer:ctx->rowOffsetBuffer offset:0 atIndex:7];
        [encoder setBuffer:ctx->rowLenBuffer offset:0 atIndex:8];
        [encoder setBuffer:ctx->resultBuffer offset:0 atIndex:9];

        /* Dispatch threads - one per row */
        const NSUInteger totalRows = (ctx->numRows > 0) ? (NSUInteger)ctx->numRows : 0;
        MTLSize gridSize = MTLSizeMake(totalRows, 1, 1);

        NSUInteger threadGroupSize = 1;
        if (pipeline != nil) {
            const NSUInteger tew = pipeline.threadExecutionWidth;
            const NSUInteger maxT = pipeline.maxTotalThreadsPerThreadgroup;
            threadGroupSize = tew * 4;
            if (threadGroupSize > maxT) threadGroupSize = maxT;
            if (threadGroupSize > 256) threadGroupSize = 256;
            if (threadGroupSize == 0) threadGroupSize = 1;
        }
        MTLSize threadgroupSize = MTLSizeMake(threadGroupSize, 1, 1);

        if (useThreadgroupKernel) {
            [encoder dispatchThreadgroups:gridSize threadsPerThreadgroup:threadgroupSize];
        } else {
            if (totalRows > 0 && threadGroupSize > totalRows) {
                threadgroupSize = MTLSizeMake(totalRows, 1, 1);
            }
            [encoder dispatchThreads:gridSize threadsPerThreadgroup:threadgroupSize];
        }
        [encoder endEncoding];

        [commandBuffer commit];
        [commandBuffer waitUntilCompleted];

        return 0;
    }
}

int metal_m2l_compute_async(MetalM2LContext* ctx) {
    if (ctx == nullptr || !ctx->setupDone) return -1;
    if (ctx->numRows <= 0) return 0;

    @autoreleasepool {
        id<MTLCommandBuffer> commandBuffer = [ctx->commandQueue commandBuffer];
        id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];
        id<MTLComputePipelineState> pipeline = nil;
        bool useThreadgroupKernel = false;

        if (ctx->threadgroupMode == 1 && ctx->pipelineFloatKahanTG != nil) {
            pipeline = ctx->pipelineFloatKahanTG;
            useThreadgroupKernel = true;
        } else if (ctx->pipelineFloatKahan != nil) {
            pipeline = ctx->pipelineFloatKahan;
        } else {
            [encoder endEncoding];
            return -1;
        }

        [encoder setComputePipelineState:pipeline];
        [encoder setBuffer:ctx->srcMmRe0Buffer offset:0 atIndex:0];
        [encoder setBuffer:ctx->srcMmIm0Buffer offset:0 atIndex:1];
        [encoder setBuffer:ctx->srcMmRe1Buffer offset:0 atIndex:2];
        [encoder setBuffer:ctx->srcMmIm1Buffer offset:0 atIndex:3];
        [encoder setBuffer:ctx->coeffReBuffer offset:0 atIndex:4];
        [encoder setBuffer:ctx->coeffImBuffer offset:0 atIndex:5];
        [encoder setBuffer:ctx->termSrcIdxBuffer offset:0 atIndex:6];
        [encoder setBuffer:ctx->rowOffsetBuffer offset:0 atIndex:7];
        [encoder setBuffer:ctx->rowLenBuffer offset:0 atIndex:8];
        [encoder setBuffer:ctx->resultBuffer offset:0 atIndex:9];

        const NSUInteger totalRows = (ctx->numRows > 0) ? (NSUInteger)ctx->numRows : 0;
        MTLSize gridSize = MTLSizeMake(totalRows, 1, 1);

        NSUInteger threadGroupSize = 1;
        if (pipeline != nil) {
            const NSUInteger tew = pipeline.threadExecutionWidth;
            const NSUInteger maxT = pipeline.maxTotalThreadsPerThreadgroup;
            threadGroupSize = tew * 4;
            if (threadGroupSize > maxT) threadGroupSize = maxT;
            if (threadGroupSize > 256) threadGroupSize = 256;
            if (threadGroupSize == 0) threadGroupSize = 1;
        }
        MTLSize threadgroupSize = MTLSizeMake(threadGroupSize, 1, 1);

        if (useThreadgroupKernel) {
            [encoder dispatchThreadgroups:gridSize threadsPerThreadgroup:threadgroupSize];
        } else {
            if (totalRows > 0 && threadGroupSize > totalRows) {
                threadgroupSize = MTLSizeMake(totalRows, 1, 1);
            }
            [encoder dispatchThreads:gridSize threadsPerThreadgroup:threadgroupSize];
        }
        [encoder endEncoding];

        [commandBuffer commit];
        ctx->pendingCommandBuffer = commandBuffer;

        return 0;
    }
}

void metal_m2l_wait(MetalM2LContext* ctx) {
    if (ctx && ctx->pendingCommandBuffer) {
        [ctx->pendingCommandBuffer waitUntilCompleted];
        ctx->pendingCommandBuffer = nil;
    }
}

void metal_m2l_destroy(MetalM2LContext* ctx) {
    if (ctx) {
        /* ARC will release the Objective-C objects when ctx is deleted */
        delete ctx;
    }
}

const char* metal_m2l_device_name(MetalM2LContext* ctx) {
    return ctx ? ctx->deviceName : "Unknown";
}

} /* extern "C" */
