---
layout: default
title: summarize
parent: IMInternalModes
grand_parent: Core
nav_order: 14
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
    factory parameters: none

  Internal-mode variables
    solved formulation: G
    diagnostic variable: F
    diagnostic relation: F_j(z) = h_j dG_j/dz(z)

  Internal-mode inner products
    F
      interior weight: 1
      endpoint terms: none
      availability: interiorOnly
      reason: For hydrostatic G modes with G=0 at both endpoints, the diagnostic F inner product is the interior F integral.
    G
      interior weight: N^2(z)/g
      endpoint terms: none
      availability: interiorOnly
      reason: The solved formulation has no endpoint metric terms, so the inner product is the interior integral only.
  ```

  For a fixed-frequency wave EVP, the diagnostic `F` inner product is
  reported as unavailable until a fixed diagnostic catalog entry is added:

  ```text
  Internal-mode inner products
    F
      interior weight: 1
      endpoint terms: none
      availability: unknown
      reason: No fixed diagnostic inner-product catalog entry is installed for this EVP and boundary combination.
  ```
