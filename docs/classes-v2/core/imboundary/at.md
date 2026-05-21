---
layout: default
title: at
parent: IMBoundary
grand_parent: Classes
nav_order: 2
mathjax: true
---

#  at

Place a boundary condition at a physical endpoint.


---

## Declaration
```matlab
 boundary = at(boundary,location,options)
```
## Parameters
+ `location`  `"surface"` or `"bottom"`
+ `options.formulation`  EVP formulation, `"F"` or `"G"`

## Returns
+ `boundary`  placed boundary condition

## Discussion

  Location-free families such as `rigid`, `free`, and `noSlip`
  resolve their variable-dependent formulas when placed. The
  assembled operator still uses coordinate $$\partial_z$$ at
  both endpoints.

  ```matlab
  condition = IMBoundary.free().at("surface", formulation="G");
  ```
