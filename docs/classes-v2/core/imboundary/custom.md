---
layout: default
title: custom
parent: IMBoundary
grand_parent: Core
nav_order: 3
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
+ `options.boundaryWeights`  endpoint weights implied by this law
+ `options.hasKnownBoundaryWeights`  true when the compatible boundary weights are known
+ `options.variable`  variable target, `"formulation"`, `"F"`, or `"G"`
+ `options.indexSign`  expected boundary-mode eigenvalue sign, `-1`, `0`, or `1`
+ `options.indexRank`  number of boundary-mode directions; currently must be `1` when `boundaryModeNumber` is supplied
+ `options.boundaryModeNumber`  explicit endpoint mode number, `-1` for surface or `-2` for bottom

## Returns
+ `boundary`  initialized boundary law

## Discussion

  The law applies `left = lambda right` at whichever endpoint
  the EVP factory places it. Both operators are written with the
  physical coordinate derivative $$\partial_z$$ at the endpoint.

  `boundaryWeights` are optional endpoint contributions associated
  with the law. Location-free weights are placed at the endpoint
  and receive the endpoint orientation sign from Green's identity.
  Explicitly located weights must match the endpoint where the law
  is resolved.
