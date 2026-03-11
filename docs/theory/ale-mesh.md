---
layout: default
title: ALE Mesh Management
---

[Japanese version (日本語)](../ja/theory/ale-mesh.html)

# ALE Mesh Management

This page describes the Arbitrary Lagrangian--Eulerian (ALE) approach used to maintain mesh quality during time-domain simulations with moving boundaries.

## Motivation

In a purely Lagrangian approach, mesh nodes on the free surface move with the fluid velocity $\nabla\phi$. While this simplifies the free-surface tracking (nodes stay on the surface by construction), it causes severe mesh distortion over time:

- Elements near a floating body stretch and compress as the body moves.
- Nodes cluster in convergence zones and spread apart in divergence zones.
- Tangled or degenerate elements lead to numerical failure.

The ALE formulation decouples the node motion from the fluid motion, allowing control over mesh quality while still tracking the free surface accurately.

## ALE Formulation

In the ALE framework, each node moves with a prescribed velocity $\mathbf{v}$ that is generally different from the fluid velocity $\nabla\phi$. The free-surface boundary conditions must be modified to account for this difference.

### Modified Material Derivative

The material derivative of any quantity $f$ following a node moving with velocity $\mathbf{v}$ (instead of the fluid velocity) becomes:

$$\frac{Df}{Dt}\bigg|_{\text{ALE}} = \frac{\partial f}{\partial t} + \mathbf{v} \cdot \nabla f.$$

More precisely, for the velocity potential on the free surface:

$$\frac{D\phi}{Dt} = \frac{\partial\phi}{\partial t} + \frac{d\boldsymbol{\chi}}{dt} \cdot \nabla\phi,$$

where $\frac{d\boldsymbol{\chi}}{dt}$ is the node velocity in the ALE frame. When $\frac{d\boldsymbol{\chi}}{dt} = \nabla\phi$, this reduces to the Lagrangian form. When $\frac{d\boldsymbol{\chi}}{dt} = \mathbf{0}$, it reduces to the Eulerian form.

### Free-Surface Conditions in ALE Form

The kinematic condition (surface tracking) becomes:

$$\frac{d\mathbf{x}}{dt} = \mathbf{v},$$

and the dynamic condition (potential evolution) becomes:

$$\frac{d\phi}{dt} = -gz - \frac{1}{2}|\nabla\phi|^2 + \mathbf{v} \cdot \nabla\phi,$$

where $\mathbf{v}$ is the chosen node velocity. The term $\mathbf{v} \cdot \nabla\phi$ is the ALE correction that accounts for the difference between node motion and fluid motion.

## Node Relocation Strategies

### Dirichlet Nodes (Free Surface)

Free-surface nodes are relocated to maintain a well-distributed mesh. The relocation strategy involves:

1. **Estimate the next surface position** using the current velocity and time step.
2. **Redistribute nodes** on the estimated surface to improve element quality, typically by Laplacian smoothing or spring-based relaxation.
3. **Compute the node velocity** $\mathbf{v}$ as the displacement divided by the time step.

The key constraint is that relocated nodes must remain on the free surface. Since the exact future surface is not known, an estimated surface (e.g., from a predictor step) is used.

### Neumann Nodes (Solid Walls)

Nodes on solid walls must remain on the wall surface. After relocation (e.g., to follow a moving body or to improve mesh quality), a correction step projects each node back onto the nearest wall surface:

1. Move the node according to the relocation algorithm.
2. **Project** the node onto the wall geometry (plane, cylinder, or general surface).
3. Adjust the node velocity to be consistent with the corrected position.

For flat walls, the projection is a simple orthogonal projection onto the plane. For curved surfaces, a nearest-point search is used.

## Mesh Quality Metric

The quality of each triangular element is monitored using the **circumradius-to-inradius ratio** (also known as the radius ratio):

$$Q = \frac{R}{r},$$

where $R$ is the circumradius and $r$ is the inradius of the triangle. For an equilateral triangle $Q = 2$, which is the optimum. As elements become elongated or degenerate, $Q$ increases without bound.

The solver monitors this metric and triggers remeshing or additional smoothing when $Q$ exceeds a threshold. This prevents the accumulation of poorly shaped elements that would degrade the accuracy of the BIE integration.

## Connection to Time Integration

The ALE node velocities are computed at each Runge--Kutta stage. The free-surface update proceeds as:

1. Solve the BVP for $\phi$ (and $\phi_t$ if a floating body is present).
2. Compute $\nabla\phi$ on the free surface.
3. Determine ALE node velocities $\mathbf{v}$.
4. Advance node positions and $\phi$ values using the ALE-modified equations.

Because the ALE correction appears only in the free-surface conditions, the [BIE formulation](boundary-integral.html) itself is unchanged.

---

**See also:** [Boundary Integral Equation](boundary-integral.html) for the BVP solved at each step, [Floating Body Dynamics](floating-body.html) for body-surface mesh considerations.
