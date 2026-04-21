---
layout: default
title: EigenmatricesForWavenumber
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 18
mathjax: true
---

#  EigenmatricesForWavenumber

Assemble the fixed-$$K$$ generalized EVP on the spectral grid.


---

## Declaration
```matlab
 [A,B] = EigenmatricesForWavenumber(self,k)
```
## Parameters
+ `self`  InternalModesSpectral instance
+ `k`  horizontal wavenumber

## Returns
+ `A`  left generalized-eigenproblem matrix
+ `B`  right generalized-eigenproblem matrix

## Discussion

              The eigenvalue equation is,
  G_{zz} - K^2 G = \frac{f_0^2 -N^2}{gh_j}G
  A = \frac{g}{f_0^2 -N^2} \left( \partial_{zz} - K^2*I \right)
  B = I
