---
layout: default
title: internalModes
parent: IMExponentialStratificationSolution
grand_parent: Analytical bases
nav_order: 6
mathjax: true
---

#  internalModes

Create an exact internal-mode basis.


---

## Declaration
```matlab
 basisSet = internalModes(solution,evp,options)
```
## Parameters
+ `evp`  internal-mode EVP
+ `options.nModes`  number of retained modes
+ `options.normalization`  active normalization
+ `options.metadata`  additional metadata

## Returns
+ `basisSet`  exact analytical internal-mode basis

## Discussion

Generalized-energy APV modes use ordinary Bessel functions on
positive-$$h$$ branches, modified Bessel functions on
negative-$$h$$ branches, and the exact integrated solution on
a zero branch. Endpoint inertia determines whether zero, one,
or two negative modes precede the zero and positive modes.
The `geostrophic` normalization remains real and positive for
every retained branch.
