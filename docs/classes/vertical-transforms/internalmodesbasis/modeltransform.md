---
layout: default
title: modelTransform
parent: InternalModesBasis
grand_parent: Classes
nav_order: 21
mathjax: true
---

#  modelTransform

Build a prefix-retained transform for model vertical grids.


---

## Declaration
```matlab
 transform = modelTransform(self,options)
```
## Parameters
+ `self`  InternalModesBasis instance
+ `options.component`  `"F"`, `"G"`, or `"both"`
+ `options.nModes`  requested retained mode count
+ `options.projectionTolerance`  leakage tolerance for prefix selection
+ `options.maxConditionNumber`  condition-number limit for prefix selection
+ `options.nTailCheck`  number of rejected tail modes checked for leakage
+ `options.nonlinearAliasingPolicy`  `"none"` or `"quadratic"`
+ `options.allowNoncanonical`  true to allow numerical wave-F projections

## Returns
+ `transform`  InternalModesTransform with prefix-retained modes

## Discussion

Projection resolvability and nonlinear de-aliasing are treated
as separate limits. If `nonlinearAliasingPolicy` is
`"quadratic"`, the transform automatically applies the
two-thirds retained-mode cap and records that diagnostic
separately from projection leakage.

$$
n_{\mathrm{retained}} =
\min(n_{\mathrm{requested}},n_{\mathrm{projection}},
n_{\mathrm{nonlinear}}).
$$
