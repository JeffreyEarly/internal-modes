---
layout: default
title: spectrum
parent: IMBasisSet
grand_parent: Core
nav_order: 29
mathjax: true
---

#  spectrum

Compute a scalar modal spectrum.


---

## Declaration
```matlab
 spectrum = spectrum(basisSet,coefficients)
```
## Parameters
+ `coefficients`  modal coefficients

## Returns
+ `spectrum`  modal spectrum

## Discussion

  For modal coefficients $$c_j$$ this returns
  $$S_j=M_{jj}|c_j|^2,$$ where $$M$$ is the full-domain scalar
  Gram matrix.
