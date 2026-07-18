---
layout: default
title: InternalModesBasis
parent: InternalModesBasis
grand_parent: Classes
nav_order: 4
mathjax: true
---

#  InternalModesBasis

Create a persisted vertical basis from canonical state.


---

## Declaration
```matlab
 basis = InternalModesBasis(options)
```
## Parameters
+ `options.z`  depth grid where modes are sampled
+ `options.F`  sampled F inverse modes
+ `options.G`  sampled G inverse modes
+ `options.h`  equivalent-depth vector
+ `options.N2`  buoyancy frequency squared sampled at `z`
+ `options.f0`  Coriolis parameter
+ `options.g`  gravitational acceleration
+ `options.D`  water-column depth
+ `options.kappa`  horizontal wavenumber metadata
+ `options.omega`  frequency metadata
+ `options.problemType`  text label for the EVP
+ `options.sourceDescription`  text label for the source solver
+ `options.componentRoleF`  F component role label
+ `options.componentRoleG`  G component role label
+ `options.forwardProjectionAvailableF`  true if F has a canonical projection
+ `options.forwardProjectionAvailableG`  true if G has a canonical projection
+ `options.orthogonalityWeightF`  F inner-product weight label
+ `options.orthogonalityWeightG`  G inner-product weight label

## Returns
+ `basis`  initialized InternalModesBasis instance

## Discussion

This constructor is intentionally a cheap state constructor for
annotated persistence. Scientific setup should usually use
`fromSolverAtFrequency`, `fromSolverAtWavenumber`, or
`fromSolvedModes`.
