---
layout: default
title: IMEigenvalueProblem
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 1
mathjax: true
---

#  IMEigenvalueProblem

Create a physical-coordinate EVP descriptor.


---

## Declaration
```matlab
 evp = IMEigenvalueProblem(options)
```
## Parameters
+ `options.name`  short EVP name
+ `options.formulation`  solved variable, `"F"` or `"G"`
+ `options.f0`  Coriolis parameter in radians per second
+ `options.g`  gravitational acceleration in meters per second squared
+ `options.leftOperator`  operator that builds the generalized-EVP matrix `A`
+ `options.rightOperator`  operator that builds the generalized-EVP matrix `B`
+ `options.innerWeights`  interior inner-product weight handles for `F` and `G`
+ `options.normalizations`  named modal normalization handles
+ `options.defaultNormalization`  natural default normalization for this EVP
+ `options.surfaceBoundary`  surface boundary law
+ `options.bottomBoundary`  bottom boundary law
+ `options.hFromEigenvalue`  equivalent-depth conversion
+ `options.hasBarotropicMode`  whether the EVP declares the barotropic mode
+ `options.indexValidationMode`  `"error"`, `"warning"`, or `"none"`
+ `options.parameters`  stored factory-specific physical inputs

## Returns
+ `evp`  initialized EVP descriptor

## Discussion

  The constructor is the low-level entry point for custom EVPs.
  Provide the solved `formulation`, the operators defining
  $$Aq=\lambda Bq,$$
  the surface and bottom boundary laws, the
  interior inner-product weights, and any normalizations that an
  `IMBasisSet` should expose. The standard factories are preferred
  for the built-in wave and hydrostatic problems because they set
  the operator, boundary, normalization, and mode-index metadata as
  one coherent contract.

  ```matlab
  left = IMOperator().plus(derivativeOrder=2);
  right = IMOperator().plus(coefficient=@(z,ctx) -ctx.N2(z)/ctx.g, derivativeOrder=0);
  evp = IMEigenvalueProblem(name="customG", formulation="G", ...
      leftOperator=left, rightOperator=right, ...
      surfaceBoundary=IMBoundary.rigid(), bottomBoundary=IMBoundary.rigid());
  ```
