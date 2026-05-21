---
layout: default
title: noSlip
parent: IMBoundary
grand_parent: Classes
nav_order: 24
mathjax: true
---

#  noSlip

Create a location-free no-slip boundary law.


---

## Declaration
```matlab
 boundary = IMBoundary.noSlip()
```
## Returns
+ `boundary`  initialized no-slip boundary condition

## Discussion

  In a `G` EVP, `noSlip` resolves to $$G_z=0$$. In an `F` EVP,
  it resolves to $$F=0$$.
