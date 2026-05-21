---
layout: default
title: rigid
parent: IMBoundary
grand_parent: Classes
nav_order: 27
mathjax: true
---

#  rigid

Create a location-free rigid boundary law.


---

## Declaration
```matlab
 boundary = IMBoundary.rigid()
```
## Returns
+ `boundary`  initialized rigid boundary condition

## Discussion

  In a `G` EVP, `rigid` resolves to $$G=0$$. In an `F` EVP, it
  resolves to $$F_z=0$$.
