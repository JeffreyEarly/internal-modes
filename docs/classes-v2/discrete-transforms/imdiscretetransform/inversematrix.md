---
layout: default
title: inverseMatrix
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 7
mathjax: true
---

#  inverseMatrix

Map retained modal coefficients back to sampled profiles.


---

## Discussion

The `inverseMatrix` is the $$n_z\times n_m$$ sampled modal basis

$$
A_{\mathrm i}=\Phi,
\qquad
(A_{\mathrm i})_{ij}=\Phi_{ij}=u_j(z_i).
$$

Row $$i$$ corresponds to `z(i)`, and column $$j$$ corresponds to
`modeNumber(j)`. The sampled modes use the normalization recorded
by `normalization`. Direct multiplication and `transformBack` are
equivalent:

```matlab
valuesByMatrix = transform.inverseMatrix*coefficients;
values = transform.transformBack(coefficients);
```
