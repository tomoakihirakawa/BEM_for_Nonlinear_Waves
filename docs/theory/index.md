---
layout: default
title: Theory
---

[Japanese version (日本語)](../ja/theory/)

# Theory

This section describes the mathematical formulation behind the BEM solver for nonlinear free-surface wave problems.

## Approach

The solver is based on **potential flow theory**: the fluid is assumed incompressible and irrotational, so the velocity field is the gradient of a scalar potential $\phi$ that satisfies the Laplace equation

$$\nabla^2 \phi = 0.$$

Rather than discretizing the entire fluid volume, we apply **Green's theorem** to convert the volume problem into a boundary integral equation (BIE) over the domain boundary $\Gamma$. This reduces the spatial dimension by one --- only the boundary surface needs to be meshed.

The boundary surface is discretized with **linear triangular elements**. At each time step a mixed boundary value problem is solved: Dirichlet conditions ($\phi$ known) on the free surface, Neumann conditions ($\partial\phi/\partial n$ known) on solid walls and the seabed.

The **Fast Multipole Method (FMM)** accelerates the BIE evaluation from $O(N^2)$ to $O(N)$, making large-scale 3D simulations feasible. The linear system is solved iteratively with GMRES.

Time integration of the free-surface position and potential uses explicit Runge--Kutta methods, with an ALE (Arbitrary Lagrangian--Eulerian) mesh update to maintain element quality.

## Contents

- [Boundary Integral Equation](boundary-integral.html) --- Green's theorem, discretization, coefficient matrices, multi-node treatment
- [Floating Body Dynamics](floating-body.html) --- 6-DOF rigid body motion, pressure integration, the $\phi_t$ problem
- [Wave Generation and Absorption](wave-generation.html) --- Piston/flap wavemakers, solitary waves, damping zones
- [ALE Mesh Management](ale-mesh.html) --- Node relocation, surface tracking, mesh quality
- [Fast Multipole Method](fmm.html) --- Spherical harmonic expansions, octree, M2L translation
