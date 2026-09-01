---
layout: default
title: rawVariable
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 22
mathjax: true
---

#  rawVariable

Evaluate raw physical `F` or `G` modes.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 values = rawVariable(basisSet,variable,z)
```
## Parameters
+ `variable`  `"F"` or `"G"`
+ `z`  physical coordinate

## Returns
+ `values`  unnormalized physical mode values

## Discussion

  This developer utility evaluates unnormalized internal-mode
  physical variables. The public `F` and `G` methods apply the
  active basis normalization after this step:
  $$F_j(z)=F_j^{\mathrm{raw}}(z)/s_j,\qquad
  G_j(z)=G_j^{\mathrm{raw}}(z)/s_j.$$
  If `variable` is the solved formulation, values come
  directly from the native modes. Otherwise the diagnostic
  variable is recovered through the EVP-owned diagnostic
  relation.
