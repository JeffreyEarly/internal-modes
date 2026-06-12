---
layout: default
title: rigid
parent: IMBoundary
grand_parent: Core
nav_order: 17
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
+ `boundary`  initialized rigid boundary law

## Discussion

  In a `G` EVP, `rigid` resolves to $$G=0$$. In an `F` EVP, it
  resolves to $$F_z=0$$.
