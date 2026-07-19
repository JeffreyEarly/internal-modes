---
layout: default
title: inverseMatrixConditionNumber
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 8
mathjax: true
---

#  inverseMatrixConditionNumber

Two-norm condition number of the sampled modal basis.


---

## Discussion

This is

$$
\kappa_2(A_{\mathrm i})
=\frac{\sigma_{\max}(A_{\mathrm i})}
{\sigma_{\min}(A_{\mathrm i})}.
$$

Large values indicate that retained mode columns are nearly
linearly dependent on the selected sample points. The definition
applies to rectangular `inverseMatrix` matrices through their
singular values.
