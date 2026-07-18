---
layout: default
title: normalizationFactors
parent: IMBasisSet
grand_parent: Core
nav_order: 17
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
The default scalar `"unity"` rule is installed on every
scalar basis set and uses
$$s_j=\sqrt{|\langle u_j,u_j\rangle|}$$ with the EVP inner
product. Custom rules are added to the basis set with
`addNormalization`, and the active rule is selected by
`basisSet.normalization`.

```matlab
factors = basisSet.normalizationFactors("unity");
uUnity = basisSet.u(z,normalization="unity");
```
