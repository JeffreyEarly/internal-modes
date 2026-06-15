---
layout: default
title: normalizationFactors
parent: IMBasisSet
grand_parent: Core
nav_order: 11
mathjax: true
---

#  normalizationFactors

Return scale factors for a normalization rule.


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

  For a requested rule name this returns the row vector
  $$s_j$$ used by evaluation methods:
  $$u_j^{\mathrm{out}}=u_j^{\mathrm{raw}}/s_j.$$
  The default scalar `"unity"` rule uses
  $$s_j=\sqrt{|\langle u_j,u_j\rangle|}$$ with the EVP inner
  product. Custom rules are created on
  `basisSet.evp.normalizationRules`, the default rule is named
  by `basisSet.evp.defaultNormalization`, and the active rule
  is selected by `basisSet.normalization`.

  ```matlab
  factors = basisSet.normalizationFactors("unity");
  uUnity = basisSet.u(z,normalization="unity");
  ```
