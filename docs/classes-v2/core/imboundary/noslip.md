---
layout: default
title: noSlip
parent: IMBoundary
grand_parent: Core
nav_order: 15
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
+ `boundary`  initialized no-slip boundary law

## Discussion

  In a `G` EVP, `noSlip` resolves to $$G_z=0$$. In an `F` EVP,
  it resolves to $$F=0$$.
