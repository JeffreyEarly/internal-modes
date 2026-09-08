---
layout: default
title: gramError
parent: IMProjection
grand_parent: Discrete transforms
nav_order: 8
mathjax: true
---

#  gramError

Magnitude-normalized active Gram operator discrepancy.


---

## Discussion

  Scaling uses the absolute diagonal of the target Gram matrix.
  A rank-deficient active sampled Gram reports infinity. A prescribed
  dual reports NaN because its system is not a component Gram.
