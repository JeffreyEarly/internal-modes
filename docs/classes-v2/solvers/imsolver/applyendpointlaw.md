---
layout: default
title: applyEndpointLaw
parent: IMSolver
grand_parent: Solvers
nav_order: 3
mathjax: true
---

#  applyEndpointLaw

Apply a resolved endpoint law to a matrix pair.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [A,B] = applyEndpointLaw(solver,A,B,endpointLaw,options)
```
## Parameters
+ `A`  left EVP matrix
+ `B`  right EVP matrix
+ `endpointLaw`  resolved endpoint law
+ `options.context`  framework coefficient context

## Returns
+ `A`  left matrix with the endpoint row applied
+ `B`  right matrix with the endpoint row applied

## Discussion

  Resolved endpoint laws replace the solver-native row associated
  with their physical location.
