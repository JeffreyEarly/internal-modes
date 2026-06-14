---
layout: default
title: negativeEigenvalueBounds
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 16
mathjax: true
---

#  negativeEigenvalueBounds

Bound negative eigenvalues using grid-level certification.


---

## Declaration
```matlab
 bounds = negativeEigenvalueBounds(evp,solver,A)
```
## Parameters
+ `solver`  canonical EVP solver
+ `A`  assembled left matrix, used for the zero-eigenvalue check

## Returns
+ `bounds`  struct with min/max counts and a reason

## Discussion
