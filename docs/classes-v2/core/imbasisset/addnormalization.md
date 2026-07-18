---
layout: default
title: addNormalization
parent: IMBasisSet
grand_parent: Core
nav_order: 2
mathjax: true
---

#  addNormalization

Add a named normalization rule.


---

## Declaration
```matlab
 basisSet = addNormalization(basisSet,name,rule)
```
## Parameters
+ `name`  normalization rule name
+ `rule`  function handle with signature `scale = rule(basisSet,iMode)`

## Returns
+ `basisSet`  basis set with the rule installed

## Discussion

`addNormalization` registers or replaces one rule on this
basis set. The rule has signature
`scale = rule(basisSet,iMode)` and returns the raw scale
factor $$s_j$$ for retained mode `iMode`. Evaluated modes use
$$u_j^{\mathrm{out}}=u_j^{\mathrm{raw}}/s_j.$$

For a constant-scaled norm,
$$s_j=C\|u_j\|_\mu,$$
use:

```matlab
C = 2;
basisSet = basisSet.addNormalization("constantScaled", @(basisSet,j) C*basisSet.innerProductNormFactor(j));
basisSet.normalization = "constantScaled";
```

For an eigenvalue-scaled norm,
$$s_j=\sqrt{|\lambda_j|}\|u_j\|_\mu,$$
use:

```matlab
basisSet = basisSet.addNormalization("eigenvalueScaled", @(basisSet,j) sqrt(abs(basisSet.eigenvalues(j)))*basisSet.innerProductNormFactor(j));
```

If `name` already exists, the rule is overwritten.
