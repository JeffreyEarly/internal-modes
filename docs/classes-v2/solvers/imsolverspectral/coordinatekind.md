---
layout: default
title: coordinateKind
parent: IMSolverSpectral
grand_parent: Solvers
nav_order: 5
mathjax: true
---

#  coordinateKind

Native coordinate kind.


---

## Discussion

  `coordinateKind` is `"z"`, `"wkb"`, or `"density"`. The `"wkb"`
stretch uses $$dx/dz=\sqrt{B}$$, where
$$B=\sqrt{(N^2+\tilde f^2-f_0^2)^2+4f_0^2\tilde f^2}$$
and $$\tilde f=2\Omega\cos\phi$$ uses the full horizontal
rotation component. At `rotationRate=0` this reduces to $$dx/dz=N$$.
This numerical monitor does not add nontraditional wave dynamics.
The `"wkb"`
  and `"density"` coordinates require an `IMInternalModes` EVP or
  `IMGeostrophicZeroAPVModes` problem because their coordinate maps
  use the problem-owned `N2` profile.
