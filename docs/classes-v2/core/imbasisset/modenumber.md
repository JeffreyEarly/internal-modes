---
layout: default
title: modeNumber
parent: IMBasisSet
grand_parent: Core
nav_order: 15
mathjax: true
---

#  modeNumber

Retained-mode labels.


---

## Discussion

  `modeNumber` is a row vector with one integer label per retained
  mode column. `modeNumber(j)` labels `eigenvalues(j)` and column
  `j` of `u(z)`, `uz(z)`, and, for internal-mode basis sets, `F(z)`
  and `G(z)`.

  Numerical solves label retained modes in the order
  $$-1,-2,\ldots,\quad 0,\quad 1,2,\ldots.$$
  Negative labels are ordinal labels for retained negative
  eigenvalues sorted by eigenvalue order; they do not identify
  surface, bottom, or any other fixed physical branch. The label
  `0` marks an inferred zero, barotropic, or null mode when one is
  retained. Positive labels mark ordinary positive/interior modes.

  ```matlab
  basisSet.modeNumber
  % [-1 0 1 2]
  ```

  See `modeSelectionDiagnostics` for the diagnostics that explain
  why negative or zero modes were retained.
