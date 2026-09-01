---
layout: default
title: certifiedDiscreteTransform
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 5
mathjax: true
---

#  certifiedDiscreteTransform

Select and independently fit a certified aligned F/G family band.


---

## Declaration
```matlab
 [transform,assessment] = certifiedDiscreteTransform(basisSet,options)
```
## Parameters
+ `options.nPoints`  exact requested mode-root point count
+ `options.z`  explicit increasing physical sample points
+ `options.variables`  requested direct channels, F and/or G
+ `options.gridDesign`  optional provenance returned by `modeRootGrid`
+ `options.gramTolerance`  per-channel normalized-Gram tolerance
+ `options.leakageTolerance`  optional rejected-mode leakage tolerance
+ `options.quadraticAliasingTolerance`  optional coupled-product tolerance
+ `options.nCheckModes`  optional rejected-mode check count

## Returns
+ `transform`  independently fitted certified aligned transform
+ `assessment`  final diagnostics, grid provenance, and count-search table

## Discussion

  Fresh family-specific weights are fitted for every count considered by
  the Gram search. This prevents a poorly resolved large candidate family
  from contaminating the assessment of a smaller MDA or APV band. Optional
  leakage and coupled quadratic policies can reduce that Gram-certified
  band; each reduction is refitted until the fitted and retained counts
  agree.
