---
layout: default
title: residualNorm
parent: IMQuadratureWeightFit
grand_parent: Discrete transforms
nav_order: 16
mathjax: true
---

#  residualNorm

Two-norm of the fitted objective residual.


---

## Discussion

This is

$$
\left\|A_{\mathrm{LS}}w-b_{\mathrm{LS}}\right\|_2.
$$

Compare it with `geometricResidualNorm` to assess improvement over
the geometric control-volume weights under the same objective.
