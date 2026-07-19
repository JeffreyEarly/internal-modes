---
layout: default
title: z
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 18
mathjax: true
---

#  z

Physical points at which profiles and modes are sampled.


---

## Discussion

`z` is an $$n_z\times1$$ column vector. Its entries define the row
ordering of `inverseMatrix`, the column ordering of
`forwardMatrix`, and the expected row ordering of values passed to
`transformForward`. Coordinates use the same units as the source
basis set, normally meters for vertical-mode problems.

```matlab
z = transform.z;
values = profile(z);
coefficients = transform.transformForward(values);
```
