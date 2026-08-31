---
layout: default
title: innerProduct
parent: IMMeanDensityAnomalyModes
grand_parent: Core
nav_order: 5
mathjax: true
---

#  innerProduct

Return the mean-density-anomaly `F` or `G` inner product.


---

## Declaration
```matlab
 spec = innerProduct(evp,variable)
```
## Parameters
+ `variable`  `"F"` or `"G"`

## Returns
+ `spec`  continuous inner-product specification

## Discussion

  `G` has the signed generalized-energy product defined by the
  canonical EVP. `F` is a diagnostic pressure shape and does
  not define a direct projection metric.
