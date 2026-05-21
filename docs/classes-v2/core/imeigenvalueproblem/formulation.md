---
layout: default
title: formulation
parent: IMEigenvalueProblem
grand_parent: Classes
nav_order: 8
mathjax: true
---

#  formulation

Solved vertical-structure formulation.


---

## Discussion

  The formulation is either `"G"` or `"F"`. The basis set solves
  this variable directly and obtains the other variable
  diagnostically through
  $$F_j=h_j\partial_zG_j$$ for `G` formulations or
  $$G_j=-gN^{-2}\partial_zF_j$$ for `F` formulations.
