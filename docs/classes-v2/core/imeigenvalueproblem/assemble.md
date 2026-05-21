---
layout: default
title: assemble
parent: IMEigenvalueProblem
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  assemble

the EVP on a solver's native basis.


---

## Declaration
```matlab
 [A,B] = assemble(evp,solver)
```
## Parameters
+ `solver`  coordinate-aware internal-mode solver

## Returns
+ `A`  left generalized-EVP matrix
+ `B`  right generalized-EVP matrix

## Discussion
