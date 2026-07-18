---
layout: default
title: internalModeAvailability
parent: IMAnalyticalSolution
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

Concrete solution families return a struct with
`isAvailable`, `reason`, `solutionKind`, `stratification`,
`supportedVariables`, `supportedInnerProducts`, and
`supportedNormalizations`.
