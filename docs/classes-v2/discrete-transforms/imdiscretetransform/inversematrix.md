---
layout: default
title: inverseMatrix
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 7
mathjax: true
---

#  inverseMatrix

Inverse transform matrix $$A_{\mathrm i}=\Phi$$.


---

## Discussion

Multiplying modal coefficients by `inverseMatrix` reconstructs
profiles on `z`. The matrix is the sampled modal basis and may be
rectangular. This is the matrix applied by `reconstruct`.
