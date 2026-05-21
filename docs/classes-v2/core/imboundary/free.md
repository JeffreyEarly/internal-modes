---
layout: default
title: free
parent: IMBoundary
grand_parent: Classes
nav_order: 12
mathjax: true
---

#  free

Create a location-free free boundary law.


---

## Declaration
```matlab
 boundary = IMBoundary.free()
```
## Returns
+ `boundary`  initialized free boundary condition

## Discussion

  In a `G` EVP, `free` resolves to $$G_z=\lambda G$$. In an `F`
  EVP, it resolves to $$F+gF_z/N^2=0$$. The derivative $$G_z$$
  or $$F_z$$ always means the coordinate derivative
  $$\partial_z$$. A placed `G` free boundary declares an endpoint
  boundary branch with mode number `-1` at the surface or `-2`
  at the bottom.
