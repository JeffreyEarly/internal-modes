---
layout: default
title: relativeGramError
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 14
mathjax: true
---

#  relativeGramError

Signed-norm-scaled relative Gram error.


---

## Discussion

With
$$S=\operatorname{diag}(|\operatorname{diag}\Gamma_0|^{-1/2})$$,
this is
$$\|S(\Gamma-\Gamma_0)S\|_2$$. When $$\Gamma_0$$ is positive
definite, this is the worst-case relative Parseval error.
