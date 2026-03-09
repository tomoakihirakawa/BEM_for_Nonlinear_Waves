---
layout: default
title: "Example: Goring (1979)"
---

# Example: Solitary Wave Generation (Goring, 1979)

This example demonstrates solitary wave generation using a piston-type wavemaker, following the method of Goring (1979).

## Reference

Goring, D.G. (1979). *Tsunamis — The propagation of long waves onto a shelf*. W.M. Keck Laboratory of Hydraulics and Water Resources, California Institute of Technology, Report No. KH-R-38.

## Setup

A 2D numerical wave tank with:
- Tank length: ~10 m
- Water depth: 0.3 m
- Piston wavemaker at the left end
- Wave absorber at the right end

The wavemaker displacement follows the solitary wave theory to generate a clean solitary wave of specified height.

## Input Files

The input files are located in `bem/input_files/` under the Goring1979 directories.

### Running

```bash
cd bem/build
./main_time_domain ../../bem/input_files/Goring1979_DT0d03_.../
```

## Expected Results

The solitary wave should propagate along the tank with:
- Wave height matching the theoretical value
- Minimal trailing oscillations
- Good agreement with wave gauge measurements at multiple positions
