---
layout: default
title: Solvers
parent: Class documentation V2
nav_order: 2
has_children: true
permalink: /classes-v2/solvers
mathjax: true
---

Reference pages for solvers that assemble canonical vertical eigenvalue problems, solve the discrete system, and return normalized basis sets.

After configuring a solver for an EVP, use `nativeQuadratureRule` to obtain its physical quadrature points and weights as one paired result:

```matlab
solver = IMSolverSpectral(nEVP=65,coordinateKind="wkb").configuredForEVP(evp);
[z,weights] = solver.nativeQuadratureRule(evp.zDomain);
```

The returned `z` values are increasing and `weights` are in the same order. For the WKB spectral coordinate, the points are the Chebyshev--Lobatto values in

$$
x(z)=\int_{z_b}^{z}\sqrt{B(s)}\,ds,\qquad B(z)=\sqrt{\left(N^2(z)+f_h^2-f_0^2\right)^2+4f_0^2f_h^2}
$$

mapped back to physical depth. Here $$f_0=2\Omega\sin\phi$$ and $$f_h=2\Omega\cos\phi$$ are derived from `rotationRate` in radians per second and `latitude` in degrees north. The bandwidth stretch reduces to the original $$N(z)$$ stretch when `rotationRate=0`. It changes the numerical coordinate while preserving the configured EVP.

The weights are the corresponding Clenshaw--Curtis value weights after the physical Jacobian $$dz/dx=1/\sqrt{B(z)}$$ is applied. The stretch must be strictly positive throughout the domain. Thus `weights.'*values` evaluates the solver's native approximation to $$\int_{z_b}^{z_s}f(z)\,dz$$ for `values = f(z)`. Keep the returned points and weights together; no additional fit is implied by this API.
