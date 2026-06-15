---
layout: default
title: partialGramMatrix
parent: IMBasisSet
grand_parent: Core
nav_order: 13
mathjax: true
---

#  partialGramMatrix

Return a partial-domain scalar Gram matrix.


---

## Declaration
```matlab
 gram = partialGramMatrix(basisSet,zMin,zMax)
```
## Parameters
+ `zMin`  lower physical bound
+ `zMax`  upper physical bound

## Returns
+ `gram`  scalar Gram matrix

## Discussion

  Interior integrals are restricted to `[zMin,zMax]`. Endpoint
  metric terms are included only when the requested interval
  contains the corresponding physical endpoint.
