---
layout: default
title: summarize
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 19
mathjax: true
---

#  summarize

Print a readable mathematical summary of this EVP.


---

## Declaration
```matlab
 summarize(evp,solver)
```
## Parameters
+ `solver`  optional solver for grid-level coefficient and mode-selection assessment

## Discussion

`summarize` prints the canonical scalar equation, physical domain,
boundary conditions, endpoint norm weights, and parameter names. When a
solver is supplied, it also samples the coefficients on the solver grid,
assembles the left matrix for the zero-mode check, and prints the
grid-level negative and zero mode assessment. It does not solve the EVP.

```matlab
evp.summarize()
```

prints output like:

```text
dirichlet

Canonical EVP
  -(d/dz)(p(z) du/dz) + q(z) u = lambda r(z) u
  z in [-1, 0]

Boundary conditions
  surface: u(surface) = 0
  bottom: u(bottom) = 0
```

```matlab
solver = IMSolverSpectral(nEVP=64);
evp.summarize(solver)
```

prints additional solver-grid output like:

```text
Solver assessment
  solver: IMSolverSpectral
  coordinate: z
  grid size: 64

Coefficient ranges on solver grid
  p: [1, 1]
  q: [0, 0]
  r: [1, 1]

Mode selection assessment
  negative modes: expected 0
  zero mode: absent
```
