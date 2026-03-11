---
layout: default
title: Floating Body Dynamics
---

[Japanese version (日本語)](../ja/theory/floating-body.html)

# Floating Body Dynamics

This page describes the coupling between the potential flow BIE solver and rigid body dynamics for floating structures.

![Boundary types with floating body](../images/bem/schematic_boundary_types_with_float.png)

*Figure: Boundary types in a wave tank with a floating body. The wetted body surface introduces additional Neumann conditions whose normal velocity depends on the body motion.*

## Rigid Body Equations of Motion

A floating body with mass $m$ and center of gravity at $\mathbf{x}_c$ obeys Newton's equations in six degrees of freedom (surge, sway, heave, roll, pitch, yaw):

$$m \frac{d\mathbf{U}_c}{dt} = \mathbf{F}_{\text{ext}} + \mathbf{F}_{\text{hydro}},$$

$$I_G \frac{d\boldsymbol{\omega}}{dt} + \boldsymbol{\omega} \times (I_G \boldsymbol{\omega}) = \mathbf{T}_{\text{ext}} + \mathbf{T}_{\text{hydro}},$$

where $\mathbf{U}_c$ is the translational velocity of the center of gravity, $\boldsymbol{\omega}$ is the angular velocity, and $I_G$ is the inertia tensor about the center of gravity. External forces $\mathbf{F}_{\text{ext}}$ may include gravity, mooring, and other loads.

### Inertia Tensor Transformation

The inertia tensor is most conveniently defined in body-fixed coordinates $I_{\text{body}}$. In the global frame it transforms as:

$$I_G = R\, I_{\text{body}}\, R^{-1},$$

where $R$ is the rotation matrix from body-fixed to global coordinates. Because $R$ is orthogonal, $R^{-1} = R^T$.

## Hydrodynamic Force and Pressure Integration

The hydrodynamic force and torque on the body are obtained by integrating the fluid pressure over the wetted surface $S_b$:

$$\mathbf{F}_{\text{hydro}} = \iint_{S_b} p\,\mathbf{n}\,dS, \qquad \mathbf{T}_{\text{hydro}} = \iint_{S_b} p\,(\mathbf{x} - \mathbf{x}_c) \times \mathbf{n}\,dS,$$

where $\mathbf{n}$ points into the fluid.

## Bernoulli Pressure

The pressure is given by the unsteady Bernoulli equation:

$$p = -\rho\left(\phi_t + \frac{1}{2}|\nabla\phi|^2 + gz\right),$$

where $\rho$ is the fluid density, $g$ is gravitational acceleration, and $\phi_t = \partial\phi/\partial t$ is the time derivative of the potential in the global (Eulerian) frame.

The terms $\nabla\phi$ and $gz$ are straightforward to evaluate once the BVP for $\phi$ is solved. The difficulty lies in computing $\phi_t$.

## The $\phi_t$ Problem

The time derivative $\phi_t$ satisfies the same Laplace equation and the same type of BIE as $\phi$:

$$\nabla^2 \phi_t = 0 \quad \text{in } \Omega,$$

$$\alpha(\mathbf{a})\,\phi_t(\mathbf{a}) = \iint_\Gamma \left[ G\,\phi_{nt} - \phi_t\,G_n \right] dS.$$

The coefficient matrices $M$ and $N$ are identical to those of the $\phi$ BIE --- only the boundary values change. This means the same FMM tree and GMRES infrastructure can be reused.

### Boundary Conditions for $\phi_t$

- **Free surface** (Dirichlet): $\phi_t$ is obtained from the dynamic free-surface condition:

$$\phi_t = -\frac{1}{2}|\nabla\phi|^2 - gz.$$

- **Solid walls and seabed** (Neumann): $\phi_{nt} = 0$ for stationary boundaries.

- **Moving body surface** (Neumann): The normal derivative $\phi_{nt}$ requires careful treatment and is the subject of the next section.

### The $\phi_{nt}$ Boundary Condition (Wu 1998)

On the wetted body surface, the kinematic condition gives $\phi_n = \mathbf{v}_b \cdot \mathbf{n}$ where $\mathbf{v}_b$ is the local velocity of the body surface. Taking the time derivative introduces the body acceleration and requires accounting for the time dependence of the surface normal. Following Wu (1998), the Neumann condition for $\phi_t$ on the body is:

$$\phi_{nt} = \dot{\mathbf{v}}_b \cdot \mathbf{n} + \text{(geometric correction terms)},$$

where $\dot{\mathbf{v}}_b$ depends on the translational acceleration $\dot{\mathbf{U}}_c$ and angular acceleration $\dot{\boldsymbol{\omega}}$ of the body. The geometric correction terms arise because the unit normal $\mathbf{n}$ itself rotates with the body.

This creates a **coupling problem**: to compute $\mathbf{F}_{\text{hydro}}$ we need $\phi_t$, which requires $\phi_{nt}$, which depends on the acceleration, which in turn depends on $\mathbf{F}_{\text{hydro}}$.

## Methods for Computing Body Acceleration

Several approaches have been proposed to resolve this implicit coupling:

### 1. Indirect Method

Compute $\phi_t$ on the body surface directly from the BIE without explicitly needing the acceleration, by rearranging the discrete system. The acceleration is then obtained from the equations of motion after the force is known.

### 2. Mode-Decomposition Method

Decompose $\phi_t$ into a part that is linear in the (unknown) acceleration and a remainder:

$$\phi_t = \sum_{k=1}^{6} \dot{q}_k\,\psi_k + \phi_t^{(0)},$$

where $\dot{q}_k$ are the six acceleration components and $\psi_k$ are unit-acceleration basis functions, each satisfying a BIE with known boundary conditions. This requires solving 7 BIEs (6 basis functions plus the particular solution), but the system matrices are identical, so only the right-hand sides differ.

### 3. Iterative Methods

**Cao iterative method:** Guess an initial acceleration, solve the $\phi_t$ BIE, compute forces, update the acceleration from the equations of motion, and iterate until convergence.

**Ma iterative method:** Similar in spirit but applies iteration at the level of the boundary condition rather than the full BIE solve.

### 4. Auxiliary Function Method (Wu and Taylor 2003)

Introduce an auxiliary potential $\psi$ that satisfies

$$\nabla^2 \psi = 0 \quad \text{in } \Omega,$$

with boundary conditions chosen so that $\psi$ absorbs the acceleration-dependent part of $\phi_{nt}$. By applying Green's theorem to $\phi_t$ and $\psi$ together, one can extract the force contribution linear in the acceleration without solving additional BIEs for each mode. This is more efficient than mode decomposition when the number of bodies is large.

### 5. Quasi-Newton Root Finding

The coupling can be viewed as a root-finding problem: find the acceleration $\mathbf{a}$ such that

$$\mathbf{a} = \frac{1}{m}\mathbf{F}_{\text{hydro}}(\mathbf{a}) + \frac{1}{m}\mathbf{F}_{\text{ext}}.$$

A quasi-Newton iteration (e.g., Broyden's method) updates an approximate Jacobian of this residual without requiring explicit derivatives. This approach is flexible and does not assume linearity in the acceleration.

## Time Integration

The coupled fluid--body system is advanced in time using explicit Runge--Kutta methods. At each stage:

1. Solve the $\phi$ BIE (mixed Dirichlet/Neumann problem).
2. Compute $\nabla\phi$ on the boundary.
3. Solve the $\phi_t$ BIE with the appropriate coupling method.
4. Evaluate hydrodynamic forces and moments.
5. Update body position, velocity, and the free-surface variables.

The [ALE mesh management](ale-mesh.html) strategy ensures that the mesh on the moving body surface remains well-conditioned throughout the simulation.

---

**See also:** [Boundary Integral Equation](boundary-integral.html) for the underlying BIE formulation, [Wave Generation](wave-generation.html) for incident wave specification.
