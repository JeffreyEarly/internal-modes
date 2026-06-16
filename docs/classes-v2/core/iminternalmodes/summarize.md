---
layout: default
title: summarize
parent: IMInternalModes
grand_parent: Core
nav_order: 15
mathjax: true
---

#  summarize

Print a readable internal-mode EVP summary.


---

## Declaration
```matlab
 summarize(evp,solver)
```
## Parameters
+ `solver`  optional solver for grid-level coefficient and mode-selection assessment

## Discussion

  `summarize` first prints the canonical EVP report, then appends the
  internal-mode interpretation: solved formulation, physical constants,
  equivalent-depth mapping, diagnostic-variable relation, factory-specific
  parameter names, and the available `F` and `G` inner products.

  ```matlab
  evp.summarize()
  ```

  prints output like:

  ```text
  Internal-mode context
    formulation: G
    f0: 0 s^-1
    g: 9.81 m s^-2
    equivalent depth: h = hFromEigenvalue(lambda)
    mode family: geostrophic
    factory parameters: none

  Internal-mode variables
    solved formulation: G
    diagnostic variable: F
    diagnostic relation: F_j(z) = h_j dG_j/dz(z)

  Internal-mode inner products
    F
      interior weight: 1
      endpoint terms: none
      inner product: known
      reason: The hydrostatic value-only catalog gives an interior-only diagnostic inner product for this boundary combination.
    G
      interior weight: N^2(z)/g
      endpoint terms: none
      inner product: known
      reason: The solved formulation has no endpoint metric terms, so the inner product is the interior integral only.

  Internal-mode normalization
    geostrophic normalization: available
      shared scale: one factor per coupled (F,G) mode
      convention: <G_j,G_j>_G = 1
      implied: <F_j,F_j>_F = h_j
  ```

  For a fixed-frequency wave EVP, the diagnostic `F` inner product is
  reported as unavailable until a wave diagnostic catalog entry is added:

  ```text
  Internal-mode inner products
    F
      interior weight: 1
      endpoint terms: none
      inner product: unavailable
      reason: The diagnostic inner product is available only for modeFamily="geostrophic" EVPs in the value-only catalog.
  ```
