---
layout: default
title: boundaryModeDescriptors
parent: IMBoundary
grand_parent: Classes
nav_order: 3
mathjax: true
---

#  boundaryModeDescriptors

Return declared endpoint boundary-mode metadata.


---

## Declaration
```matlab
 descriptors = boundaryModeDescriptors(boundary)
```
## Returns
+ `descriptors`  structure array with `modeNumber` and `indexSign`

## Discussion

  Boundary-mode numbers are endpoint labels used for retained
  modes. They are selected using `indexSign`, but the sign of the
  eigenvalue and the physical mode number are intentionally
  separate concepts.
