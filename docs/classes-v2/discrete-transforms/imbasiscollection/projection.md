---
layout: default
title: projection
parent: IMBasisCollection
grand_parent: Discrete transforms
nav_order: 9
mathjax: true
---

#  projection

Construct a fixed projection for one requested page and variable.


---

## Declaration
```matlab
 [projection,recipe] = projection(collection,z,weights,options)
```
## Parameters
+ `z`  physical sample column
+ `weights`  fixed quadrature weights aligned with z
+ `options.page`  requested-page index, not distinct-solve index
+ `options.variable`  scalar u or aligned F/G variable
+ `options.columns`  explicit column indices; empty means all

## Returns
+ `projection`  IMProjection with the requested columns
+ `recipe`  basis-owned metric, capabilities, and provenance

## Discussion

  This operation preserves the supplied samples, weights, and explicit
  columns. It performs no quadrature fitting or acceptance policy. The
  continuous basis supplies the signed pairing and positive error metric.
