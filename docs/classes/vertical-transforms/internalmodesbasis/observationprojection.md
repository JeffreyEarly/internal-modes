---
layout: default
title: observationProjection
parent: InternalModesBasis
grand_parent: Classes
nav_order: 23
mathjax: true
---

#  observationProjection

Build a vertical projection for an arbitrary observation grid.


---

## Declaration
```matlab
 projection = observationProjection(self,observation,options)
```
## Parameters
+ `self`  InternalModesBasis instance
+ `observation`  observation depths or explicit observation matrix
+ `options.component`  `"F"` or `"G"`
+ `options.nModes`  number of candidate modes
+ `options.weights`  observation weights
+ `options.rankTolerance`  relative QR pivot tolerance
+ `options.maxConditionNumber`  retained Gram condition-number limit

## Returns
+ `projection`  InternalModesProjection for the observation operator

## Discussion

The observation problem uses the sampled matrix
$$B=H\Phi$$, where $$H$$ is either a point-sampling grid or an
explicit observation operator and $$\Phi$$ is the inverse mode
matrix. Rank-revealing QR selects a resolvable, possibly
non-contiguous set of modes.
