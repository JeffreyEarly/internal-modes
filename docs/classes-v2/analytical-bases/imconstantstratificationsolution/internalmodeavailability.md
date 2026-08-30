---
layout: default
title: internalModeAvailability
parent: IMConstantStratificationSolution
grand_parent: Analytical bases
nav_order: 4
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

For a `"geostrophicAPVModes"` descriptor, availability also
verifies constant stratification, the canonical APV
coefficients, both endpoint conditions, and the direct APV
metadata contract. Qualifying hydrostatic families report
the volume-only `depth` normalization.
