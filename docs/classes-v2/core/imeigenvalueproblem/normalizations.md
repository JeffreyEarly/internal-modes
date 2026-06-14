---
layout: default
title: normalizations
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 17
mathjax: true
---

#  normalizations

Named normalization rules.


---

## Discussion

  Each field stores a function handle with signature
  `scale = rule(basisSet,iMode)`. The returned value is the raw
  scale $$s_j$$ for one mode, and basis-set evaluation divides all
  variables for that mode by $$s_j$$. The `unity` rule is supplied
  automatically when omitted. Internal-mode factories add rules for
  `geostrophic`, `kConstant`, `omegaConstant`, `wMax`, `uMax`, and
  `surfacePressure`.

  ```matlab
  normalizations.unity = @(basisSet,iMode) ...
      basisSet.innerProductNormFactor(iMode);
  evp = IMEigenvalueProblem(normalizations=normalizations, ...
      defaultNormalization=Normalization.unity);
  ```
