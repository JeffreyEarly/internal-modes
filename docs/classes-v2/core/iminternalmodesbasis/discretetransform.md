---
layout: default
title: discreteTransform
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 8
mathjax: true
---

#  discreteTransform

Build an aligned internal-mode F/G transform through the compatibility API.


---

## Declaration
```matlab
 [transform,assessment] = discreteTransform(basisSet,options)
```
## Parameters
+ `options.nPoints`  exact requested mode-root point count
+ `options.z`  explicit increasing sample locations
+ `options.weights`  optional explicit shared weights
+ `options.nModes`  optional strict family prefix for explicit z
+ `options.variables`  requested direct channels, F and/or G
+ `options.gramTolerance`  per-channel normalized-Gram tolerance
+ `options.leakageTolerance`  optional rejected-mode tolerance
+ `options.quadraticAliasingTolerance`  optional coupled-product tolerance
+ `options.nCheckModes`  optional rejected-mode check count

## Returns
+ `transform`  retained aligned transform
+ `assessment`  family diagnostics and, for certified construction, grid and search provenance

## Discussion

  The point rule is shared, but each requested variable retains its own
  sampled metric, active columns, target Gram matrix, and forward
  projection. Omitted `variables` selects every directly representable
  channel in canonical order `F`, `G`. Point-limited construction still
  uses roots of the next mode in the EVP's solved formulation.

  With neither `weights` nor `nModes`, this delegates to
  `certifiedDiscreteTransform`, which independently refits family weights
  while selecting the retained count. Prefer that named method in new code.
  `fitDiscreteTransform(z=z,modeCount=N)` is the strict exact-band API, and
  `modeRootGrid` makes the source of shared APV/MDA points explicit.

  Supplying weights or the legacy `nModes` name retains fixed-rule prefix
  assessment. The Gram policy must pass independently for every requested
  channel. Optional leakage uses same-variable rejected modes. Coupled
  quadratic aliasing assesses `FF->F` and `GG->F` when F is enabled, and
  `FG->G` when G is enabled.
