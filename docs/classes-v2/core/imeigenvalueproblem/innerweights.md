---
layout: default
title: innerWeights
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 17
mathjax: true
---

#  innerWeights

Inner-product weights for `F` and `G`.


---

## Discussion

  `innerWeights.F` and `innerWeights.G` are function handles with
  signature `w = weight(z,ctx)`. Each handle returns the interior
  weight in
  $$\langle X_i,X_j\rangle_w=\int w(z)X_i(z)X_j(z)\,dz.$$
  Endpoint contributions belong to the boundary laws so that each
  normalization can include the same endpoint convention as the EVP.
