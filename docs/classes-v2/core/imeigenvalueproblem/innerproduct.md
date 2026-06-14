---
layout: default
title: innerProduct
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 11
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
  endpoint metric terms implied by active endpoint conditions.
