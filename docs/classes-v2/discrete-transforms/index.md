---
layout: default
title: Discrete transforms
parent: Class documentation V2
nav_order: 4
has_children: true
permalink: /classes-v2/discrete-transforms
mathjax: true
---

Reference pages for constructing scalar Galerkin transforms, fitting quadrature increments, transforming sampled profiles forward to modal coefficients, transforming coefficients back to sample space, and assessing discrete Parseval accuracy.

The workflow begins with a solved `IMBasisSet`: `quadraturePoints` proposes a mode-root grid, `fitQuadrature` compares fitted and geometric increments, and `IMDiscreteTransform` stores `forwardMatrix` and `inverseMatrix`. Use `transformForward` to compute retained modal coefficients and `transformBack` to return those coefficients to the transform sample points.
