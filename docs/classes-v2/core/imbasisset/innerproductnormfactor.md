---
layout: default
title: innerProductNormFactor
parent: IMBasisSet
grand_parent: Core
nav_order: 8
mathjax: true
---

#  innerProductNormFactor

Return the scalar inner-product norm factor.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 factor = innerProductNormFactor(basisSet,iMode)
```
## Parameters
+ `iMode`  retained mode index

## Returns
+ `factor`  raw scalar inner-product scale factor

## Discussion

  This developer utility returns the raw scale factor used by
  the default `"unity"` normalization rule before any modal
  normalization has been applied:
  $$s_j=\sqrt{|\langle u_j,u_j\rangle_\mu|}.$$
  The scalar inner product includes the interior weight and
  any canonical endpoint terms prepared by `endpointGramTerms`.
