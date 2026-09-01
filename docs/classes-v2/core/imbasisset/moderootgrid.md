---
layout: default
title: modeRootGrid
parent: IMBasisSet
grand_parent: Core
nav_order: 16
mathjax: true
---

#  modeRootGrid

Design a mode-root grid and report how it was generated.


---

## Declaration
```matlab
 [z,gridDesign] = modeRootGrid(basisSet,options)
```
## Parameters
+ `options.nPoints`  exact requested physical point count
+ `options.modeCount`  leading mode count represented by the generating mode

## Returns
+ `z`  increasing physical grid
+ `gridDesign`  scalar provenance struct

## Discussion

  Specify either an exact physical point count or the represented leading
  mode count. The returned `gridDesign` records the source EVP and family,
  generating variable and mode, and the interpretation of the points for
  internal-mode `G` structures.

  ```matlab
  [z,gridDesign] = basisSet.modeRootGrid(nPoints=128);
  gridDesign.generatingVariable
  gridDesign.interpretationForG
  ```
