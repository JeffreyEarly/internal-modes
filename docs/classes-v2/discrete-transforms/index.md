---
layout: default
title: Discrete transforms
parent: Class documentation V2
nav_order: 4
has_children: true
permalink: /classes-v2/discrete-transforms
mathjax: true
---

Reference pages for constructing mode-root point grids, finding quadrature weights for fixed points, constructing scalar Galerkin transforms, transforming sampled profiles forward to modal coefficients, transforming coefficients back to sample space, and assessing discrete Parseval accuracy.

The workflow begins with a solved `IMBasisSet`: `pointsFromModeRoots` returns both endpoints and the interior roots of the next selected mode. Given those fixed points, `discreteTransform` fits the quadrature weights, constructs the scalar transform, and optionally returns the associated `IMQuadratureWeightFit`. `IMDiscreteTransform` stores `forwardMatrix` and `inverseMatrix`; use `transformForward` to compute retained modal coefficients and `transformBack` to return those coefficients to the transform sample points.

```matlab
z = basisSet.pointsFromModeRoots(nModes=8);
[transform,weightFit] = basisSet.discreteTransform(z=z,nModes=8);
```

Use `quadratureWeightsForPoints` directly when the fitted weights or a custom linear objective are the primary concern:

```matlab
[weights,weightFit] = basisSet.quadratureWeightsForPoints(z=z,nModes=8);
transform = basisSet.discreteTransform(z=z,weights=weights,nModes=8);
```
