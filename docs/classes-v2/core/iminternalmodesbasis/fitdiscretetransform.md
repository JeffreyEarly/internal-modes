---
layout: default
title: fitDiscreteTransform
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 10
mathjax: true
---

#  fitDiscreteTransform

Fit and certify one exact aligned F/G family band.


---

## Declaration
```matlab
 [transform,assessment] = fitDiscreteTransform(basisSet,options)
```
## Parameters
+ `options.z`  increasing physical sample points
+ `options.modeCount`  exact number of leading family modes
+ `options.variables`  requested direct channels, F and/or G
+ `options.gridDesign`  optional provenance returned by `modeRootGrid`
+ `options.gramTolerance`  per-channel normalized-Gram tolerance
+ `options.leakageTolerance`  optional rejected-mode leakage tolerance
+ `options.quadraticAliasingTolerance`  optional coupled-product tolerance
+ `options.nCheckModes`  optional rejected-mode check count

## Returns
+ `transform`  exact fitted aligned transform
+ `assessment`  exact-band diagnostics and grid provenance

## Discussion

  `modeCount` is explicit and strict: the weights are fitted to exactly
  that many aligned family columns, and all requested variables and enabled
  policies must accept the complete band. Use `certifiedDiscreteTransform`
  to select a retained count automatically.
