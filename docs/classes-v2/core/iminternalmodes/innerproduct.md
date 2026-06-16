---
layout: default
title: innerProduct
parent: IMInternalModes
grand_parent: Core
nav_order: 12
mathjax: true
---

#  innerProduct

Return the `F` or `G` inner-product recipe.


---

## Declaration
```matlab
 spec = innerProduct(evp,variable)
```
## Parameters
+ `variable`  optional variable name, `"F"` or `"G"`

## Returns
+ `spec`  struct with interior and endpoint inner-product terms

## Discussion

  For `G`, the interior weight is $$N^2/g$$. For `F`, the
  interior weight is one. The returned struct has fields
  `variable`, `interiorWeight`, `surfaceWeights`,
  `bottomWeights`, `endpointInnerProductTerms`,
  `hasInnerProduct`, and `reason`. `hasInnerProduct` is true
  when the variable has a known inner product. When it is false, Gram
  matrices, spectra, and inner-product normalization for that
  variable throw `IMInternalModesBasis:UnavailableInnerProduct`.
  Diagnostic variables use the value-only hydrostatic endpoint
  catalog only when `modeFamily` is `"geostrophic"` and a
  catalog row is known; other diagnostic inner products are
  unavailable until a family catalog is added. Endpoint
  inner-product terms from the catalog have the form
  $$\alpha_\ell F_i(z_\ell)F_j(z_\ell)$$ or
  $$\alpha_\ell G_i(z_\ell)G_j(z_\ell),$$
  where $$z_\ell$$ is the bottom or surface endpoint. The
  variable used in the endpoint term is stored as
  `term.variable`, so a `G` inner product can contain an
  endpoint term involving `F`, and conversely.
