---
layout: default
title: normalizationFactors
parent: IMBasisSet
grand_parent: Core
nav_order: 13
mathjax: true
---

#  normalizationFactors

Return factors for a normalization convention.


---

## Declaration
```matlab
 factors = normalizationFactors(basisSet,normalization)
```
## Parameters
+ `normalization`  normalization convention

## Returns
+ `factors`  row vector of positive scale factors

## Discussion

  For a requested convention this returns the row vector
  $$s_j$$ used by evaluation methods:
  $$u_j^{\mathrm{out}}=u_j^{\mathrm{raw}}/s_j.$$
  The default scalar `unity` convention uses
  $$s_j=\sqrt{|\langle u_j,u_j\rangle|}$$ with the EVP inner
  product.

  ```matlab
  factors = basisSet.normalizationFactors(Normalization.unity);
  uUnity = basisSet.u(z,normalization=Normalization.unity);
  ```
