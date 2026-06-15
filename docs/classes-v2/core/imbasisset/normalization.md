---
layout: default
title: normalization
parent: IMBasisSet
grand_parent: Core
nav_order: 10
mathjax: true
---

#  normalization

Active normalization rule name.


---

## Discussion

  This string selects a field in `basisSet.evp.normalizationRules`;
  create custom rules on the EVP, not on the basis set.
  The selected rule returns the scale factor $$s_j$$ used by `u`,
  `uz`, and Gram-matrix methods. Passing `normalization=...` to an
  evaluation method overrides this property for that call.
