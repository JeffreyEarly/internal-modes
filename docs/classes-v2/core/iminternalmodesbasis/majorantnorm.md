---
layout: default
title: majorantNorm
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 17
mathjax: true
---

#  majorantNorm

Return the positive Hilbert-majorant norm of modal coefficients.


---

## Declaration
```matlab
 value = majorantNorm(basisSet,coefficients,options)
```
## Parameters
+ `coefficients`  one coefficient per retained mode
+ `options.variable`  `"F"` or `"G"`
+ `options.zBounds`  integration bounds `[zMin zMax]`

## Returns
+ `value`  positive scalar norm

## Discussion

  For coefficient vector $$c$$ this returns
  $$\sqrt{c^*M_+c}$$, where $$M_+$$ is
  `majorantGramMatrix`. Do not replace this with
  $$\sqrt{|c^*Mc|}$$ for a signed Gram matrix $$M$$; the latter
  can vanish for a nonzero state and is not a norm.
