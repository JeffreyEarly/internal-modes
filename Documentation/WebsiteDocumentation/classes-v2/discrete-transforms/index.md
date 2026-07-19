---
layout: default
title: Discrete transforms
parent: Class documentation V2
nav_order: 4
has_children: true
permalink: /classes-v2/discrete-transforms
mathjax: true
---

Reference pages for choosing quadrature points, finding quadrature weights for fixed points, constructing scalar Galerkin transforms, transforming sampled profiles forward to modal coefficients, transforming coefficients back to sample space, and assessing discrete Parseval accuracy.

The workflow begins with a solved `IMBasisSet`: `quadraturePoints` proposes a mode-root grid, and `quadratureWeightsForPoints` finds weights for any fixed grid. `IMQuadratureWeightFit` optionally compares those fitted algebraic weights with geometric control-volume weights. `IMDiscreteTransform` stores `forwardMatrix` and `inverseMatrix`; use `transformForward` to compute retained modal coefficients and `transformBack` to return those coefficients to the transform sample points.

```matlab
z = basisSet.quadraturePoints(nModes=8);
weights = basisSet.quadratureWeightsForPoints(z=z,nModes=8);
transform = basisSet.discreteTransform(z=z,weights=weights,nModes=8);
```
