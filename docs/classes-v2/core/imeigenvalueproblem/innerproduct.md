---
layout: default
title: innerProduct
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 12
mathjax: true
---

#  innerProduct

Return the scalar inner-product recipe.


---

## Declaration
```matlab
 spec = innerProduct(evp,variable)
```
## Parameters
+ `variable`  scalar variable name; only `"u"` is accepted

## Returns
+ `spec`  struct with interior and endpoint metric terms

## Discussion

  The canonical basis set uses `r` in the interior and the
  endpoint metric terms implied by active endpoint conditions:
  $$M_{ij}=\int r u_i u_j\,dz+
  \sum_\ell \gamma_\ell L_\ell[u_i]L_\ell[u_j].$$
  The returned struct has fields `variable`, `interiorWeight`,
  `surfaceWeights`, `bottomWeights`, and
  `hasKnownBoundaryWeights`.
