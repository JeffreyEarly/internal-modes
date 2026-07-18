---
layout: default
title: canonicalBoundary
parent: IMHydrostaticBoundaryCondition
grand_parent: Core
nav_order: 5
mathjax: true
---

#  canonicalBoundary

Convert to a canonical scalar boundary condition.


---

## Declaration
```matlab
 boundary = canonicalBoundary(law,formulation="G",g=9.81)
```
## Parameters
+ `options.formulation`  target solved variable, `"F"` or `"G"`
+ `options.g`  gravitational acceleration

## Returns
+ `boundary`  canonical boundary condition

## Discussion

`canonicalBoundary` converts the physical hydrostatic law to
the `IMBoundaryCondition` object required by a canonical
scalar EVP. Use `formulation="F"` for hydrostatic `F` EVPs
and `formulation="G"` for hydrostatic `G` EVPs. The
conversion follows the constructor-shaped rules in the class
overview.
