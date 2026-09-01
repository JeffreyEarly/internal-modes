---
layout: default
title: fitDiscreteTransform
parent: IMBasisSet
grand_parent: Core
nav_order: 9
mathjax: true
---

#  fitDiscreteTransform

Fit and certify one exact modal band on caller-selected points.


---

## Declaration
```matlab
 [transform,assessment] = fitDiscreteTransform(basisSet,options)
```
## Parameters
+ `options.z`  increasing physical sample points
+ `options.modeCount`  exact number of leading modes to fit and certify
+ `options.gridDesign`  optional provenance returned by `modeRootGrid`
+ `options.gramTolerance`  normalized-Gram operator-error tolerance
+ `options.leakageTolerance`  optional rejected-mode leakage tolerance
+ `options.quadraticAliasingTolerance`  optional quadratic-aliasing tolerance
+ `options.nCheckModes`  optional rejected-mode check count

## Returns
+ `transform`  exact fitted transform
+ `assessment`  exact-band diagnostics and grid provenance

## Discussion

  This is the strict, diagnostic API. `modeCount` is the number of leading
  family modes whose quadrature weights are fitted. The same exact band is
  then required to pass every enabled policy. Construction throws rather
  than silently returning a shorter prefix.

  Use `certifiedDiscreteTransform` when the retained count should be chosen
  automatically. Use `modeRootGrid` when the physical points should be
  designed from a particular modal family before fitting one or more
  families independently on that shared grid.
