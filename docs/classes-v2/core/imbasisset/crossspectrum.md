---
layout: default
title: crossSpectrum
parent: IMBasisSet
grand_parent: Core
nav_order: 2
mathjax: true
---

#  crossSpectrum

Compute a scalar modal cross-spectrum.


---

## Declaration
```matlab
 spectrum = crossSpectrum(basisSet,coefficientsA,coefficientsB)
```
## Parameters
+ `coefficientsA`  first modal coefficients
+ `coefficientsB`  second modal coefficients

## Returns
+ `spectrum`  modal cross-spectrum

## Discussion

  For modal coefficient vectors $$a_j$$ and $$b_j$$ this
  returns $$S_j=M_{jj}\operatorname{Re}(a_j b_j^*)$$.
