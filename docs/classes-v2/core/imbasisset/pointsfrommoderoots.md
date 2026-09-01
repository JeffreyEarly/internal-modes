---
layout: default
title: pointsFromModeRoots
parent: IMBasisSet
grand_parent: Core
nav_order: 24
mathjax: true
---

#  pointsFromModeRoots

Return physical endpoints and roots of the next selected mode.


---

## Declaration
```matlab
 z = pointsFromModeRoots(basisSet,options)
```
## Parameters
+ `options.nModes`  number of leading selected basis columns represented by the grid

## Returns
+ `z`  increasing physical grid containing both endpoints and the generating mode's interior roots

## Discussion

  Let `nModes=N` select the first $$N$$ columns of the ordered basis. The
  generating mode is the next selected column,

  $$
  u_{\mathrm{gen}}(z)=u_{\mathrm{selected},N+1}(z).
  $$

  This is a column index in the selected basis, not a physical mode label.
  In particular, `modeNumber(N+1)` need not equal $$N+1$$ when negative or
  zero modes precede the positive modes. The returned physical grid is

  $$
  \mathbf{z}_{\mathrm{root}}=\{z_b\}\cup
  \left\{z\in(z_b,z_s):u_{\mathrm{gen}}(z)=0\right\}\cup\{z_s\}.
  $$

  The interior roots are converted to physical $$z$$ coordinates, sorted,
  and deduplicated. Both physical endpoints are always included. If column
  $$N+1$$ is not already retained, the stored solver and EVP obtain one
  auxiliary mode and verify that the new solve reproduces the original
  selected prefix; the basis set itself is not changed.

  The result is a deterministic mode-root candidate grid, not a complete
  quadrature rule. Use `quadratureWeightsForPoints` to fit weights and assess
  how accurately the resulting rule reproduces the retained modal Gram
  matrix.

  ```matlab
  z = basisSet.pointsFromModeRoots(nModes=8);
  weights = basisSet.quadratureWeightsForPoints(z=z,nModes=8);
  transform = basisSet.discreteTransform(z=z,weights=weights,nModes=8);
  ```
