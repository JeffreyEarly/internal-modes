---
layout: default
title: internalModeAvailability
parent: IMExponentialStratificationSolution
grand_parent: Analytical bases
nav_order: 5
mathjax: true
---

#  internalModeAvailability

Report whether exact internal modes are available.


---

## Declaration
```matlab
 availability = internalModeAvailability(solution,evp)
```
## Parameters
+ `evp`  internal-mode EVP

## Returns
+ `availability`  availability report struct

## Discussion

For a `"geostrophicAPVModes"` EVP, availability additionally
verifies exponential stratification, the canonical APV
coefficients, both endpoint conditions, and the `g0`, `gd`,
and `surfaceBoundary` parameters.
