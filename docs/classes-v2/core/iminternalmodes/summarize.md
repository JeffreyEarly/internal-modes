---
layout: default
title: summarize
parent: IMInternalModes
grand_parent: Core
nav_order: 12
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
  equivalent-depth mapping, and factory-specific parameter names.

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
  ```
