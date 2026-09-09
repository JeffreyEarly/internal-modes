---
layout: default
title: rotateWithPencil
parent: IMGeostrophicZeroAPVModesBasis
grand_parent: Core
nav_order: 19
mathjax: true
---

#  rotateWithPencil

Apply a custom symmetric matrix-pencil rotation.


---

## Declaration
```matlab
 basisSet = rotateWithPencil(boundaryModes,options)
```
## Parameters
+ `options.name`  custom rotation name
+ `options.leftMatrix`  symmetric left matrix pages
+ `options.rightMatrix`  symmetric right matrix pages

## Returns
+ `basisSet`  custom-pencil rotated basis

## Discussion

  `leftMatrix` and `rightMatrix` are expressed in canonical
  boundary coordinates. Two-dimensional matrices are broadcast
  across wavenumber pages; otherwise their third dimension must
  equal `nK`.
  The pagewise pencil is

  $$
  \mathsf L\mathbf c^a=\lambda_a\mathsf R\mathbf c^a,
  \qquad
  |(\mathbf c^a)^T\mathsf R\mathbf c^a|=1.
  $$

  ```matlab
  Hg = boundaryModes.generalizedEnergyMatrix(g0=-0.035,gd=0.01);
  customModes = boundaryModes.rotateWithPencil(name="custom",leftMatrix=Hg,rightMatrix=boundaryModes.endpointResponseMetric);
  ```
