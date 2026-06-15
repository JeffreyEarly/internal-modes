---
layout: default
title: defaultNormalization
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 5
mathjax: true
---

#  defaultNormalization

Default normalization rule name.


---

## Discussion

  If a basis set is created without an explicit normalization,
  `defaultNormalization` becomes the active `basisSet.normalization`.
  The value is a string naming a field in `normalizationRules`.
  Evaluated modes are always raw modes divided by a per-mode scale,
  $$u_j^{\mathrm{out}}(z)=u_j^{\mathrm{raw}}(z)/s_j.$$
