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
quantities before subtracting `objectiveTarget`.
