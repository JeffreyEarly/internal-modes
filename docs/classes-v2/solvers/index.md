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
x(z)=\int_{z_b}^{z}N(s)\,ds
$$

mapped back to physical depth. The weights are the corresponding Clenshaw--Curtis value weights after the physical Jacobian $$dz/dx=1/N(z)$$ is applied. Thus `weights.'*values` evaluates the solver's native approximation to $$\int_{z_b}^{z_s}f(z)\,dz$$ for `values = f(z)`. Keep the returned points and weights together; no additional fit is implied by this API.
