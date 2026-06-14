---
layout: default
title: defaultNormalization
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 6
mathjax: true
---

#  defaultNormalization

Natural default normalization.


---

## Discussion

  If a basis set is created without an explicit normalization,
  `defaultNormalization` becomes the active `basisSet.normalization`.
  Evaluated modes are always raw modes divided by a per-mode scale,
  $$u_j^{\mathrm{out}}(z)=u_j^{\mathrm{raw}}(z)/s_j.$$
