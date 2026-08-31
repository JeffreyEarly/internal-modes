---
layout: default
title: Analytical bases
parent: Class documentation V2
nav_order: 3
has_children: true
permalink: /classes-v2/analytical-bases
mathjax: true
---

Reference pages for constant- and exponential-stratification exact solution families, internal-mode bases, and geostrophic zero-APV boundary modes.

## Capability follows the concrete class

`IMAnalyticalSolution` stores the common domain, Coriolis parameter, gravity, and stratification interface. It does not advertise optional solution operations. A concrete analytical family exposes a formula by implementing its construction method directly:

```matlab
solution = IMConstantStratificationSolution(N0=5.2e-3,zDomain=[-4000 0],f0=1e-4);
internalModes = solution.internalModes(evp,nModes=8);
boundaryModes = solution.geostrophicZeroAPVModesAtWavenumber(1e-4);
```

Invalid EVPs and unsupported boundary requests are rejected by the concrete construction method. There is no separate availability-report step.

## Exact generalized-energy APV modes

Both analytical families support the public `IMInternalModes.geostrophicAPVModes` descriptor. They apply either the free-surface or rigid-lid convention, handle signed finite, zero (Dirichlet), and positive-infinite endpoint accelerations, and return modes in ascending eigenvalue order $$1/h$$: up to two negative modes, an optional exact zero mode with `h=Inf`, then positive modes.

For $$N^2=N_0^2$$, the exact catalog uses hyperbolic negative branches, an affine zero branch, and trigonometric positive branches. For $$N^2(z)=N_0^2\exp(2z/b)$$, it uses the corresponding exact modified-Bessel, integrated-zero, and ordinary-Bessel branches. APV bases use the volume-only `depth` normalization by default. Roots, normalized determinant residuals, branch classifications, endpoint inertia, and the direct `g0`, `gd`, and `surfaceBoundary` contract are exposed in basis metadata and mode-selection diagnostics.

## Exact geostrophic zero-APV modes

For each positive horizontal wavenumber $$k$$, the zero-APV structure satisfies

$$
\frac{\partial}{\partial z}\left(\frac{f_0^2}{N^2}\frac{\partial F}{\partial z}\right)-k^2F=0,
\qquad
G=-\frac{g}{N^2}\frac{\partial F}{\partial z}.
$$

The canonical response coordinates are

$$
B_s[F]=
\begin{cases}
G(0)-F(0),&\text{free surface},\\
G(0),&\text{rigid lid},
\end{cases}
\qquad
B_d[F]=G(z_b).
$$

They are related to the physical endpoint displacements by $$\eta_s=(f_0/g)B_s$$ and $$\eta_d=(f_0/g)B_d$$. A two-endpoint exact basis is normalized by

$$
\mathbf B[F_0^{\mathrm{sur}}]=(1,0)^T,
\qquad
\mathbf B[F_0^{\mathrm{bot}}]=(0,1)^T.
$$

Surface-only and bottom-only requests retain the corresponding column and enforce zero response at the omitted endpoint:

```matlab
surfaceModes = solution.geostrophicZeroAPVModesAtWavenumber(k,endpoints="surface",surfaceBoundary="freeSurface");
bottomModes = solution.geostrophicZeroAPVModesAtWavenumber(k,endpoints="bottom",surfaceBoundary="rigidLid");
bothModes = solution.geostrophicZeroAPVModesAtWavenumber(k,endpoints=["surface" "bottom"]);
```

`IMAnalyticalGeostrophicZeroAPVModesBasis` evaluates exact `F(z)` and `G(z)` arrays without numerical differentiation or a numerical solver. Its endpoint-response, physical-energy, surface-buoyancy, bottom-buoyancy, generalized-energy, and rotation APIs match the numerical `IMGeostrophicZeroAPVModesBasis` contract.

The traditional one-boundary rigid-lid profiles are members of this broader exact solution space. The V2 API names the physical zero-APV family and its canonical endpoint coordinates rather than treating those profiles as a separate basis abstraction.
