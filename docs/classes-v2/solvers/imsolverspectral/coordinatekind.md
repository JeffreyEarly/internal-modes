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
  and `"density"` coordinates require an `IMInternalModes` EVP
  because their coordinate maps use the EVP-owned `N2` profile.
