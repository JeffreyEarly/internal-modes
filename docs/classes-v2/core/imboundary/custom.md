---
layout: default
title: custom
parent: IMBoundary
grand_parent: Classes
nav_order: 6
mathjax: true
---

#  custom

Create a custom location-free operator boundary law.


---

## Declaration
```matlab
 boundary = IMBoundary.custom(options)
```
## Parameters
+ `options.left`  left-side boundary functional
+ `options.right`  right-side eigenvalue functional
+ `options.innerProductTerms`  boundary inner-product terms
+ `options.hasKnownInnerProductTerms`  true when the compatible boundary inner-product terms are known
+ `options.variable`  variable target, `"formulation"`, `"F"`, or `"G"`
+ `options.indexSign`  expected boundary-mode eigenvalue sign, `-1`, `0`, or `1`
+ `options.indexRank`  number of boundary-mode directions; currently must be `1` when `boundaryModeNumber` is supplied
+ `options.boundaryModeNumber`  explicit endpoint mode number, `-1` for surface or `-2` for bottom

## Returns
+ `boundary`  initialized boundary condition

## Discussion

  The law applies `left = lambda right` at whichever endpoint
  the EVP factory places it. Both operators are written with the
  physical coordinate derivative $$\partial_z$$ at the endpoint.

  `innerProductTerms` are optional boundary trace products
  associated with the law. If a term is passed with an empty
  location, it is placed at the resolved endpoint and receives
  the endpoint orientation sign from Green's identity. Explicitly
  located terms are treated as final and are not reoriented.
