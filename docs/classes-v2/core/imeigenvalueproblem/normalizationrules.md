---
layout: default
title: normalizationRules
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 15
mathjax: true
---

#  normalizationRules

Named normalization rule table.


---

## Discussion

  `normalizationRules` is where custom canonical normalization
  rules are created. Each field stores a function handle with
  signature
  `scale = rule(basisSet,iMode)`. The returned value is the raw
  scale $$s_j$$ for one mode, and basis-set evaluation divides all
  variables for that mode by $$s_j$$. `basisSet.normalization`
  selects a field from this table, and
  `basisSet.normalizationFactors(...)` evaluates that rule for all
  retained modes. The canonical `unity` rule is supplied
  automatically when omitted:
  $$s_j=\|u_j\|_\mu.$$

  A rule can depend on constants captured by the function handle,
  on the basis set, or on mode-specific quantities such as
  `basisSet.eigenvalues(iMode)`. For example, a constant-scaled
  norm can use
  $$s_j=C\|u_j\|_\mu,$$
  while an eigenvalue-scaled norm can use
  $$s_j=\sqrt{|\lambda_j|}\|u_j\|_\mu.$$

  ```matlab
  C = 2;
  rules.constantScaled = @(basisSet,iMode) ...
      C*basisSet.innerProductNormFactor(iMode);
  rules.eigenvalueScaled = @(basisSet,iMode) ...
      sqrt(abs(basisSet.eigenvalues(iMode))) * ...
      basisSet.innerProductNormFactor(iMode);
  evp = IMEigenvalueProblem(normalizationRules=rules, ...
      defaultNormalization="constantScaled");
  basisSet = solver.solveEVP(evp,nModes=4);
  basisSet.normalization = "eigenvalueScaled";
  factors = basisSet.normalizationFactors();
  ```

  The internal-mode `Normalization` enum is a convenience for
  `IMInternalModes` conventions. Generic canonical EVPs use string
  rule names.
