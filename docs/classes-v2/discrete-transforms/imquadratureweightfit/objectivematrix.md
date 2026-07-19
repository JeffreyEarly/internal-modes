---
layout: default
title: objectiveMatrix
parent: IMQuadratureWeightFit
grand_parent: Discrete transforms
nav_order: 12
mathjax: true
---

#  objectiveMatrix

Least-squares matrix $$A_{\mathrm{LS}}$$.

> Developer documentation: this item describes internal implementation details.


---

## Discussion

It has one column per fixed point. Multiplying this matrix by a
physical quadrature-weight vector produces the modeled objective
quantities before subtracting `objectiveTarget`. For the default
normalized Gram objective with $$n_m$$ retained modes, it has
$$n_m(n_m+1)/2$$ rows ordered as
$$(1,1),(1,2),\ldots,(1,n_m),(2,2),\ldots,(n_m,n_m)$$. Its
off-diagonal rows already include the $$\sqrt{2}$$ factor needed
for the residual norm to equal the full normalized Gram Frobenius
norm. Custom objectives may use any finite row count.
