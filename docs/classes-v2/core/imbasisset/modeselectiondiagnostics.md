---
layout: default
title: modeSelectionDiagnostics
parent: IMBasisSet
grand_parent: Core
nav_order: 9
mathjax: true
---

#  modeSelectionDiagnostics

Mode-selection diagnostics.


---

## Discussion

  This is the diagnostics struct returned by
  `evp.modeSelectionDiagnostics` when the solver selected and
  labeled retained modes. Numerical solves use it to record
  negative-mode bounds and zero-mode status.
