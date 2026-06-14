---
layout: default
title: endpointNumeratorMatrix
parent: IMBoundaryCondition
grand_parent: Core
nav_order: 8
mathjax: true
---

#  endpointNumeratorMatrix

Return the active-endpoint numerator matrix.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 H = endpointNumeratorMatrix(boundary,location)
```
## Discussion

  The vector is `[u; p*u_z]`.
