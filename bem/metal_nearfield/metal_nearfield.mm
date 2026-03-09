/*
 * metal_nearfield.mm
 * Objective-C++ implementation of Metal-accelerated nearfield direct integration
 *
 * Compile with Apple Clang:
 *   clang++ -c metal_nearfield.mm -o metal_nearfield.o -std=c++17 -fobjc-arc
 *   clang++ -shared -o libmetal_nearfield.dylib metal_nearfield.o -framework Metal -framework Foundation
 */

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#include "metal_nearfield.h"
#include <cstring>

struct MetalNearfieldContext {
    id<MTLDevice> device;
    id<MTLCommandQueue> commandQueue;
    id<MTLComputePipelineState> pipeline;
    id<MTLLibrary> library;

    /* Geometry buffers (set during setup, updated for ALE) */
    id<MTLBuffer> faceDataBuffer;
    id<MTLBuffer> targetPosBuffer;

    /* Structure buffers (set during setup, static) */
    id<MTLBuffer> groupFaceOffsetBuffer;
    id<MTLBuffer> groupDofCountBuffer;
    id<MTLBuffer> targetGroupIdBuffer;
    id<MTLBuffer> targetVtxIdBuffer;
    id<MTLBuffer> targetDofOffsetBuffer;

    /* Gauss quadrature constant buffers */
    id<MTLBuffer> gaussNBuffer;
    id<MTLBuffer> gaussWBuffer;

    /* Scalar parameters */
    id<MTLBuffer> nearRegionSqBuffer;
    id<MTLBuffer> numTargetsBuffer;

    /* Result buffers (Kahan pairs) */
    id<MTLBuffer> resultPhiHiBuffer;
    id<MTLBuffer> resultPhiLoBuffer;
    id<MTLBuffer> resultPhinHiBuffer;
    id<MTLBuffer> resultPhinLoBuffer;

    int32_t numTargets;
    int32_t numFaces;
    int32_t numCellGroups;
    int32_t totalResultSlots;
    bool setupDone;

    char deviceName[256];
};

extern "C" {

int metal_nearfield_is_available(void) {
    @autoreleasepool {
        id<MTLDevice> device = MTLCreateSystemDefaultDevice();
        return (device != nil) ? 1 : 0;
    }
}

MetalNearfieldContext* metal_nearfield_init(void) {
    @autoreleasepool {
        auto* ctx = new MetalNearfieldContext();
        memset(ctx, 0, sizeof(MetalNearfieldContext));

        ctx->device = MTLCreateSystemDefaultDevice();
        if (ctx->device == nil) {
            delete ctx;
            return nullptr;
        }

        const char* name = [[ctx->device name] UTF8String];
        strncpy(ctx->deviceName, name, sizeof(ctx->deviceName) - 1);

        ctx->commandQueue = [ctx->device newCommandQueue];
        if (ctx->commandQueue == nil) {
            delete ctx;
            return nullptr;
        }

        /* Load Metal shader library */
        NSError* error = nil;

        /* Try precompiled metallib first */
        NSString* libPath = @"metal_nearfield.metallib";
        NSURL* libURL = [NSURL fileURLWithPath:libPath];
        ctx->library = [ctx->device newLibraryWithURL:libURL error:&error];

        if (ctx->library == nil) {
            /* Fallback: compile from source */
            NSFileManager* fm = [NSFileManager defaultManager];
            NSArray* searchPaths = @[
                @".",
                @"./metal_nearfield",
                [[NSBundle mainBundle] bundlePath],
                @"/Users/tomoaki/Library/CloudStorage/Dropbox/code/cpp/builds/build_bem/metal_nearfield"
            ];

            NSString* foundPath = nil;
            for (NSString* dir in searchPaths) {
                NSString* fullPath = [dir stringByAppendingPathComponent:@"metal_nearfield.metal"];
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
                    options.mathMode = MTLMathModeSafe;
                    ctx->library = [ctx->device newLibraryWithSource:source
                                                             options:options
                                                               error:&error];
                }
            }
        }

        if (ctx->library == nil) {
            NSLog(@"[Metal NF] Failed to load library: %@", error);
            delete ctx;
            return nullptr;
        }

        /* Create pipeline */
        id<MTLFunction> func = [ctx->library newFunctionWithName:@"nearfield_linear_kahan"];
        if (func == nil) {
            NSLog(@"[Metal NF] Kernel function not found");
            delete ctx;
            return nullptr;
        }

        ctx->pipeline = [ctx->device newComputePipelineStateWithFunction:func error:&error];
        if (ctx->pipeline == nil) {
            NSLog(@"[Metal NF] Failed to create pipeline: %@", error);
            delete ctx;
            return nullptr;
        }

        ctx->setupDone = false;
        NSLog(@"[Metal NF] Initialized on: %s", ctx->deviceName);
        return ctx;
    }
}

int metal_nearfield_setup(
    MetalNearfieldContext* ctx,
    int32_t num_targets,
    int32_t num_faces,
    int32_t num_cell_groups,
    const MetalFaceData* face_data,
    const int32_t* group_face_offset,
    const int32_t* group_dof_count,
    const float* target_pos,
    const int32_t* target_group_id,
    const int32_t* target_vtx_id,
    const int32_t* target_dof_offset,
    int32_t total_result_slots,
    float near_region_sq)
{
    if (ctx == nullptr || ctx->pipeline == nil) return -1;

    @autoreleasepool {
        ctx->numTargets = num_targets;
        ctx->numFaces = num_faces;
        ctx->numCellGroups = num_cell_groups;
        ctx->totalResultSlots = total_result_slots;

        /* Create geometry buffers */
        ctx->faceDataBuffer = [ctx->device newBufferWithBytes:face_data
                                                       length:num_faces * sizeof(MetalFaceData)
                                                      options:MTLResourceStorageModeShared];
        ctx->targetPosBuffer = [ctx->device newBufferWithBytes:target_pos
                                                        length:num_targets * 3 * sizeof(float)
                                                       options:MTLResourceStorageModeShared];

        /* Create structure buffers */
        ctx->groupFaceOffsetBuffer = [ctx->device newBufferWithBytes:group_face_offset
                                                              length:(num_cell_groups + 1) * sizeof(int32_t)
                                                             options:MTLResourceStorageModeShared];
        ctx->groupDofCountBuffer = [ctx->device newBufferWithBytes:group_dof_count
                                                            length:num_cell_groups * sizeof(int32_t)
                                                           options:MTLResourceStorageModeShared];
        ctx->targetGroupIdBuffer = [ctx->device newBufferWithBytes:target_group_id
                                                            length:num_targets * sizeof(int32_t)
                                                           options:MTLResourceStorageModeShared];
        ctx->targetVtxIdBuffer = [ctx->device newBufferWithBytes:target_vtx_id
                                                          length:num_targets * sizeof(int32_t)
                                                         options:MTLResourceStorageModeShared];
        ctx->targetDofOffsetBuffer = [ctx->device newBufferWithBytes:target_dof_offset
                                                              length:num_targets * sizeof(int32_t)
                                                             options:MTLResourceStorageModeShared];

        /* Create Gauss quadrature constant buffers */
        /* 25-point ModTriShape<3> table - same values as GaussTableSIMD::gw5x5 */
        {
            static const double gw5[5] = {0.046910077030668, 0.230765344947158, 0.5,
                                           0.769234655052841, 0.953089922969332};
            static const double ww5[5] = {0.118463442528095, 0.239314335249683, 0.284444444444444,
                                          0.239314335249683, 0.118463442528095};
            float N[75]; /* 25 * 3 */
            float W[25];
            for (int i = 0; i < 5; i++) {
                for (int j = 0; j < 5; j++) {
                    double t0 = gw5[i], t1 = gw5[j];
                    double ww = ww5[i] * ww5[j];
                    int k = i * 5 + j;
                    /* ModTriShape<3>(t0, t1) = {t0, t1*(1-t0), (t0-1)*(t1-1)} */
                    N[k*3+0] = (float)t0;
                    N[k*3+1] = (float)(t1 * (1.0 - t0));
                    N[k*3+2] = (float)((t0 - 1.0) * (t1 - 1.0));
                    W[k] = (float)(ww * (1.0 - t0));
                }
            }
            ctx->gaussNBuffer = [ctx->device newBufferWithBytes:N
                                                         length:75 * sizeof(float)
                                                        options:MTLResourceStorageModeShared];
            ctx->gaussWBuffer = [ctx->device newBufferWithBytes:W
                                                         length:25 * sizeof(float)
                                                        options:MTLResourceStorageModeShared];
        }

        /* Scalar parameter buffers */
        ctx->nearRegionSqBuffer = [ctx->device newBufferWithBytes:&near_region_sq
                                                           length:sizeof(float)
                                                          options:MTLResourceStorageModeShared];
        ctx->numTargetsBuffer = [ctx->device newBufferWithBytes:&num_targets
                                                        length:sizeof(int32_t)
                                                       options:MTLResourceStorageModeShared];

        /* Result buffers (Kahan pairs) */
        size_t resultSize = total_result_slots * sizeof(float);
        ctx->resultPhiHiBuffer = [ctx->device newBufferWithLength:resultSize
                                                          options:MTLResourceStorageModeShared];
        ctx->resultPhiLoBuffer = [ctx->device newBufferWithLength:resultSize
                                                          options:MTLResourceStorageModeShared];
        ctx->resultPhinHiBuffer = [ctx->device newBufferWithLength:resultSize
                                                           options:MTLResourceStorageModeShared];
        ctx->resultPhinLoBuffer = [ctx->device newBufferWithLength:resultSize
                                                           options:MTLResourceStorageModeShared];

        ctx->setupDone = true;

        NSLog(@"[Metal NF] Setup: %d targets, %d faces, %d groups, %d result slots (%.1f MB)",
              num_targets, num_faces, num_cell_groups, total_result_slots,
              (double)(total_result_slots * 4 * sizeof(float)) / (1024.0 * 1024.0));
        return 0;
    }
}

int metal_nearfield_update_faces(MetalNearfieldContext* ctx, const MetalFaceData* face_data) {
    if (ctx == nullptr || !ctx->setupDone) return -1;
    memcpy([ctx->faceDataBuffer contents], face_data, ctx->numFaces * sizeof(MetalFaceData));
    return 0;
}

int metal_nearfield_update_targets(MetalNearfieldContext* ctx, const float* target_pos) {
    if (ctx == nullptr || !ctx->setupDone) return -1;
    memcpy([ctx->targetPosBuffer contents], target_pos, ctx->numTargets * 3 * sizeof(float));
    return 0;
}

int metal_nearfield_compute(MetalNearfieldContext* ctx) {
    if (ctx == nullptr || !ctx->setupDone) return -1;
    if (ctx->numTargets <= 0) return 0;

    @autoreleasepool {
        id<MTLCommandBuffer> commandBuffer = [ctx->commandQueue commandBuffer];
        id<MTLComputeCommandEncoder> encoder = [commandBuffer computeCommandEncoder];

        [encoder setComputePipelineState:ctx->pipeline];
        [encoder setBuffer:ctx->faceDataBuffer        offset:0 atIndex:0];
        [encoder setBuffer:ctx->groupFaceOffsetBuffer  offset:0 atIndex:1];
        [encoder setBuffer:ctx->groupDofCountBuffer    offset:0 atIndex:2];
        [encoder setBuffer:ctx->targetPosBuffer        offset:0 atIndex:3];
        [encoder setBuffer:ctx->targetGroupIdBuffer    offset:0 atIndex:4];
        [encoder setBuffer:ctx->targetVtxIdBuffer      offset:0 atIndex:5];
        [encoder setBuffer:ctx->targetDofOffsetBuffer  offset:0 atIndex:6];
        [encoder setBuffer:ctx->resultPhiHiBuffer      offset:0 atIndex:7];
        [encoder setBuffer:ctx->resultPhiLoBuffer      offset:0 atIndex:8];
        [encoder setBuffer:ctx->resultPhinHiBuffer     offset:0 atIndex:9];
        [encoder setBuffer:ctx->resultPhinLoBuffer     offset:0 atIndex:10];
        [encoder setBuffer:ctx->gaussNBuffer           offset:0 atIndex:11];
        [encoder setBuffer:ctx->gaussWBuffer           offset:0 atIndex:12];
        [encoder setBuffer:ctx->nearRegionSqBuffer     offset:0 atIndex:13];
        [encoder setBuffer:ctx->numTargetsBuffer       offset:0 atIndex:14];

        NSUInteger totalThreads = (NSUInteger)ctx->numTargets;
        NSUInteger tew = ctx->pipeline.threadExecutionWidth;
        NSUInteger tgSize = tew;
        if (tgSize > ctx->pipeline.maxTotalThreadsPerThreadgroup)
            tgSize = ctx->pipeline.maxTotalThreadsPerThreadgroup;
        if (tgSize > totalThreads)
            tgSize = totalThreads;
        if (tgSize == 0)
            tgSize = 1;

        MTLSize gridSize = MTLSizeMake(totalThreads, 1, 1);
        MTLSize threadgroupSize = MTLSizeMake(tgSize, 1, 1);
        [encoder dispatchThreads:gridSize threadsPerThreadgroup:threadgroupSize];

        [encoder endEncoding];
        [commandBuffer commit];
        [commandBuffer waitUntilCompleted];

        return 0;
    }
}

float* metal_nearfield_get_phi_hi(MetalNearfieldContext* ctx) {
    return ctx && ctx->resultPhiHiBuffer ? (float*)[ctx->resultPhiHiBuffer contents] : nullptr;
}

float* metal_nearfield_get_phi_lo(MetalNearfieldContext* ctx) {
    return ctx && ctx->resultPhiLoBuffer ? (float*)[ctx->resultPhiLoBuffer contents] : nullptr;
}

float* metal_nearfield_get_phin_hi(MetalNearfieldContext* ctx) {
    return ctx && ctx->resultPhinHiBuffer ? (float*)[ctx->resultPhinHiBuffer contents] : nullptr;
}

float* metal_nearfield_get_phin_lo(MetalNearfieldContext* ctx) {
    return ctx && ctx->resultPhinLoBuffer ? (float*)[ctx->resultPhinLoBuffer contents] : nullptr;
}

void metal_nearfield_destroy(MetalNearfieldContext* ctx) {
    if (ctx) {
        delete ctx;
    }
}

const char* metal_nearfield_device_name(MetalNearfieldContext* ctx) {
    return ctx ? ctx->deviceName : "Unknown";
}

} /* extern "C" */
