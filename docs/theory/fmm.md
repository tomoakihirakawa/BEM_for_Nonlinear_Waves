---
layout: default
title: Fast Multipole Method
---

[Japanese version (日本語)](../ja/theory/fmm.html)

# Fast Multipole Method

This page describes the Fast Multipole Method (FMM) used to accelerate the boundary integral equation (BIE) evaluation.

## Motivation

The [BIE discretization](boundary-integral.html) produces a dense $N \times N$ linear system, where $N$ is the number of boundary nodes. Direct solution costs $O(N^3)$ with LU factorization, and even a single matrix-vector product costs $O(N^2)$. For the 3D wave simulations targeted by this solver, $N$ can reach tens of thousands or more, making direct methods impractical.

The FMM reduces the cost of a matrix-vector product to $O(N)$ by exploiting the rapid decay of interactions between well-separated groups of source and evaluation points. Combined with the GMRES iterative solver, this gives an overall solution cost of $O(N)$ per iteration.

## Spherical Harmonic Expansions

The Green's function $G(\mathbf{x}, \mathbf{a}) = -1/(4\pi|\mathbf{x} - \mathbf{a}|)$ can be expanded in spherical harmonics when $\mathbf{x}$ and $\mathbf{a}$ are referred to a common center. The expansion separates the dependence on source and evaluation point:

- **Multipole expansion (M):** Groups the contributions from nearby sources into a compact representation centered at the group's center. Valid far from the sources.
- **Local expansion (L):** Converts the far-field influence into a Taylor-like series centered near the evaluation points. Valid close to the center.

The truncation order $p$ of the expansion controls the accuracy. Higher $p$ gives better accuracy at increased cost per translation.

## Adaptive Octree Decomposition

The FMM organizes the boundary nodes into a hierarchical spatial tree:

1. **Root cell** encloses the entire domain.
2. Each cell is recursively subdivided into eight children (octree) until the number of nodes per leaf cell falls below a threshold.
3. The tree depth adapts to the local node density --- regions with many nodes are refined more.

The tree structure defines "near" and "far" interactions:

- **Near field:** Cells at the same tree level that are adjacent. These interactions are computed directly (no approximation).
- **Far field:** Cells that are well-separated. These interactions are computed via the multipole/local expansions.

The separation criterion (MAC --- multipole acceptance criterion) ensures that the expansion error is controlled.

## Translation Operators

The FMM uses three translation operators to propagate information through the tree:

### Particle to Multipole (P2M)
Converts source values (charges, dipoles from the BIE integrands) into multipole expansion coefficients at the center of each leaf cell.

### Multipole to Multipole (M2M)
Shifts a multipole expansion from a child cell's center to its parent cell's center. Applied bottom-up through the tree.

### Multipole to Local (M2L)
Converts a multipole expansion at one cell's center into a local expansion at a well-separated cell's center. This is the most expensive operation and where different algorithmic variants differ.

### Local to Local (L2L)
Shifts a local expansion from a parent cell's center to its child cell's center. Applied top-down through the tree.

### Local to Particle (L2P)
Evaluates the local expansion at each node within a leaf cell to obtain the far-field contribution to that node.

## M2L Translation Methods

The solver provides several implementations of the M2L operator, selected at compile time:

| Method | Description |
|--------|-------------|
| `SimpleM2L` | Direct spherical harmonic translation. Straightforward but $O(p^4)$ per translation. Reliable baseline. |
| `FourierM2L_double` | Fourier-based rotation--translation--rotation decomposition. Reduces cost to $O(p^3)$ per translation using double precision. |
| `FourierM2L_DoubleDouble` | Same algorithm as above but uses double-double arithmetic for improved accuracy at high expansion orders. |
| `PlaneWaveM2L` | Plane-wave (exponential) expansion method. Asymptotically the fastest at $O(p^2)$ per translation but requires careful implementation of the plane-wave quadrature. |

For most simulations, `SimpleM2L` or `FourierM2L_double` provides a good balance of accuracy and performance.

## Connection to GMRES

The FMM provides a fast method to compute the matrix-vector product $A\mathbf{x}$, but does not produce the matrix $A$ explicitly. This rules out direct solvers and requires an iterative method.

GMRES (Generalized Minimal Residual) is used because:

- It works with the matrix only through matrix-vector products.
- It minimizes the residual norm over the Krylov subspace at each iteration.
- It handles the non-symmetric systems that arise from the BIE formulation.

A typical solve requires 10--30 GMRES iterations for the $\phi$ BIE with a tolerance of $10^{-9}$. The $\phi_t$ BIE (see [Floating Body Dynamics](floating-body.html)) uses the same FMM tree and typically converges in a similar number of iterations.

The total cost per time step is therefore $O(kN)$ where $k$ is the number of GMRES iterations --- effectively linear in $N$ since $k$ is usually independent of $N$.

---

**See also:** [Boundary Integral Equation](boundary-integral.html) for the linear system being solved, [Floating Body Dynamics](floating-body.html) for the $\phi_t$ BIE that reuses the same FMM infrastructure.
