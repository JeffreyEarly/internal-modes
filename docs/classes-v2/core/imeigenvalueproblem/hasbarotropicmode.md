---
layout: default
title: hasBarotropicMode
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 12
mathjax: true
---

#  hasBarotropicMode

Whether the EVP declares the barotropic mode.


---

## Discussion

  When `true`, the mode-index policy expects the depth-uniform
  zero-eigenvalue barotropic mode
  $$F_0(z)=1,\qquad G_0(z)=0.$$
  The barotropic mode is selected after boundary-index modes and
  before positive baroclinic modes. Hydrostatic `F` modes declare this
  mode; wave-mode and hydrostatic `G` EVPs do not.
