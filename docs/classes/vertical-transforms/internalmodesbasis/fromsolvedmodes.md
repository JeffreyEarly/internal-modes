---
layout: default
title: fromSolvedModes
parent: InternalModesBasis
grand_parent: Classes
nav_order: 13
mathjax: true
---

#  fromSolvedModes

Create a basis from already-solved mode arrays.


---

## Declaration
```matlab
 basis = InternalModesBasis.fromSolvedModes(z,F,G,h,options)
```
## Parameters
+ `z`  depth grid
+ `F`  sampled F inverse modes
+ `G`  sampled G inverse modes
+ `h`  equivalent-depth vector
+ `options.N2`  buoyancy frequency squared sampled at `z`
+ `options.f0`  Coriolis parameter
+ `options.g`  gravitational acceleration
+ `options.kappa`  horizontal wavenumber metadata
+ `options.omega`  frequency metadata
+ `options.problemType`  text label for the EVP

## Returns
+ `basis`  InternalModesBasis containing the supplied modes

## Discussion

Use this factory when another workflow has already computed
sampled inverse modes and equivalent depths.
