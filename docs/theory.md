---
layout: default
title: Theory
---

[🇯🇵 日本語](ja/theory.html)

# Theoretical Background

## Potential Flow Formulation

The fluid is assumed to be inviscid, incompressible, and irrotational. The velocity field is described by a velocity potential $\phi$ satisfying Laplace's equation:

$$\nabla^2 \phi = 0$$

## Boundary Integral Equation

Using Green's second identity, the boundary integral equation is:

$$c(\mathbf{x})\phi(\mathbf{x}) = \int_S \left[ G(\mathbf{x}, \mathbf{y}) \frac{\partial \phi}{\partial n}(\mathbf{y}) - \phi(\mathbf{y}) \frac{\partial G}{\partial n}(\mathbf{x}, \mathbf{y}) \right] dS(\mathbf{y})$$

where $G(\mathbf{x}, \mathbf{y}) = \frac{1}{4\pi|\mathbf{x} - \mathbf{y}|}$ is the free-space Green's function and $c(\mathbf{x})$ is the solid angle coefficient.

## Mixed Eulerian-Lagrangian (MEL) Method

The free surface is tracked in a Lagrangian frame. At each time step:

1. **Solve BVP**: Given $\phi$ on the free surface and $\partial\phi/\partial n$ on solid boundaries, solve the boundary integral equation for the unknown values.
2. **Compute velocities**: Calculate $\nabla\phi$ on the free surface.
3. **Update free surface**: Advance the free-surface position and potential using the kinematic and dynamic boundary conditions:

$$\frac{D\mathbf{x}}{Dt} = \nabla\phi, \quad \frac{D\phi}{Dt} = -gz + \frac{1}{2}|\nabla\phi|^2$$

## Fast Multipole Method (FMM)

The FMM accelerates the evaluation of boundary integrals from O(N²) to O(N). This implementation uses:

- Spherical harmonic expansions for multipole and local coefficients
- Multipole-to-Local (M2L) translation operators
- Adaptive octree spatial decomposition

### Available M2L Methods

| Method | Description |
|--------|-------------|
| `SimpleM2L` | Direct spherical harmonic translation (default) |
| `FourierM2L_double` | Fourier-based rotation with double precision |
| `FourierM2L_DoubleDouble` | Fourier-based with double-double precision |
| `PlaneWaveM2L` | Plane-wave expansion method |

## Arbitrary Lagrangian-Eulerian (ALE) Mesh Motion

To prevent mesh distortion, an ALE approach is used:

1. Compute Lagrangian node velocities
2. Apply Laplacian smoothing to redistribute nodes
3. Correct the material derivative for the ALE velocity

## Wave Generation

Supported wave generation methods:

- **Piston/flap wavemaker** with prescribed motion
- **Solitary wave** generation (Goring, 1979)
- **Regular waves** via linear/second-order wavemaker theory
- **Absorbing boundaries** using numerical beaches
