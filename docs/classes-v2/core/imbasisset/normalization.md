---
layout: default
title: normalization
parent: IMBasisSet
grand_parent: Core
nav_order: 11
mathjax: true
---

#  normalization

Active normalization rule name.


---

## Discussion

  This string selects a rule in the basis-set normalization
  registry. Create custom rules with `addNormalization`.
  The selected rule returns the scale factor $$s_j$$ used by `u`,
  `uz`, and Gram-matrix methods. Passing `normalization=name` to an
  evaluation method overrides this property for that call.
