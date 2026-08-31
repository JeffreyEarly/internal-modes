---
layout: default
title: transformSourceForward
parent: IMGeostrophicTransform
grand_parent: Discrete transforms
nav_order: 9
mathjax: true
---

#  transformSourceForward

Project generic vorticity and displacement sources.


---

## Declaration
```matlab
 [APVSourceCoefficients,zeroAPVSourceCoefficients] = transformSourceForward(transform,options)
```
## Parameters
+ `options.vorticitySource`  sampled vorticity-source pages
+ `options.displacementSource`  sampled displacement-source pages
+ `options.zeroAPVCoordinates`  output zero-APV coordinates

## Returns
+ `APVSourceCoefficients`  APV source-coefficient pages
+ `zeroAPVSourceCoefficients`  zero-APV source-coefficient pages

## Discussion

  This method uses the APV transform's
  `modeProjectionFunctional` operation rather than
  `transformForward`. The source equations require raw modal
  pairings before the sampled Gram solve.

  $$
  S_q^j=\frac{1}{D}\int F_jS_\omega\,dz
  -\frac{f_0}{D}\mathcal G_q^j[S_\eta].
  $$

  Zero-APV source coefficients satisfy

  $$
  \widehat{\mathsf H}_g^\kappa\mathbf S_0^\kappa
  =\mathbf p_0^\kappa,\qquad
  \widehat{\mathsf H}_g^\kappa=\frac{2\kappa^2}{D}\mathsf H_g^\kappa.
  $$
