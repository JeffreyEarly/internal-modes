---
layout: default
title: certifiedDiscreteTransform
parent: IMBasisSet
grand_parent: Core
nav_order: 3
mathjax: true
---

#  certifiedDiscreteTransform

Select and independently fit a certified scalar modal band.


---

## Declaration
```matlab
 [transform,assessment] = certifiedDiscreteTransform(basisSet,options)
```
## Parameters
+ `options.nPoints`  exact requested mode-root point count
+ `options.z`  explicit increasing physical sample points
+ `options.gridDesign`  optional provenance returned by `modeRootGrid`
+ `options.gramTolerance`  normalized-Gram operator-error tolerance
+ `options.leakageTolerance`  optional rejected-mode leakage tolerance
+ `options.quadraticAliasingTolerance`  optional quadratic-aliasing tolerance
+ `options.nCheckModes`  optional rejected-mode check count

## Returns
+ `transform`  independently fitted certified transform
+ `assessment`  final diagnostics, grid provenance, and count-search table

## Discussion

  Each candidate count receives a fresh quadrature-weight fit to exactly
  that band. The search therefore never judges a short prefix using weights
  fitted to a larger, poorly resolved family. The largest Gram-certified
  count is selected first. Enabled leakage or quadratic policies may then
  reduce the count; every reduction is refitted and reassessed until the
  fitted and retained counts agree.

  `nPoints` asks this basis to design a mode-root grid. Explicit `z` may be
  accompanied by the `gridDesign` returned by another family's
  `modeRootGrid`, making shared-grid provenance inspectable.
