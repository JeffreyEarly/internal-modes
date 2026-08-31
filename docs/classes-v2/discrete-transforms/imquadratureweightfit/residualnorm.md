---
layout: default
title: residualNorm
parent: IMQuadratureWeightFit
grand_parent: Discrete transforms
nav_order: 16
mathjax: true
---

#  residualNorm

Norm of the fitted objective residual.


---

## Discussion

  This is

  $$
  \left\|A_{\mathrm{LS}}w-b_{\mathrm{LS}}\right\|_2.
  $$

  For `objectiveName="normalizedGramFrobenius"`, the least-squares
  residual is the vectorization of the normalized Gram mismatch, so

  $$
  \texttt{residualNorm}=\|E(w)\|_{\mathrm F}.
  $$

  For a custom objective, this property retains the generic meaning
  $$\|Aw-b\|_2$$. Compare it with `geometricResidualNorm` only under
  the same objective.
