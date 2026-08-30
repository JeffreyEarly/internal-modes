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
and `"density"` coordinates require an `IMInternalModes` EVP or
`IMGeostrophicZeroAPVModes` problem because their coordinate maps
use the problem-owned `N2` profile.
