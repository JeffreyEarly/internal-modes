---
layout: default
title: IMEigenvalueProblem
parent: IMEigenvalueProblem
grand_parent: Classes
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
+ `options.f0`  Coriolis parameter
+ `options.g`  gravitational acceleration
+ `options.leftOperator`  left physical operator
+ `options.rightOperator`  right physical operator
+ `options.innerWeights`  interior inner-product weights for `F` and `G`
+ `options.normalizations`  named modal normalization rules
+ `options.defaultNormalization`  natural default normalization for this EVP
+ `options.boundaryConditions`  placed boundary conditions
+ `options.hFromEigenvalue`  equivalent-depth conversion
+ `options.nNullModes`  number of true null modes
+ `options.indexValidationMode`  `"error"`, `"warning"`, or `"none"`
+ `options.parameters`  stored factory-specific physical inputs

## Returns
+ `evp`  initialized EVP descriptor

## Discussion
