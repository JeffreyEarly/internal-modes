---
layout: default
title: G
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 2
mathjax: true
---

#  G

Evaluate `G` modes.


---

## Declaration
```matlab
 G = G(basisSet,z,options)
```
## Parameters
+ `z`  physical coordinate
+ `options.normalization`  normalization rule name or enum value

## Returns
+ `G`  evaluated `G` modes

## Discussion

If the EVP formulation is `G`, this evaluates the solved
canonical variable. If the formulation is `F`, `G` is recovered
by `evp.GfromFz`. The default relation is
$$G_j(z)=-\frac{g}{N^2(z)}
\frac{\partial F_j}{\partial z}(z),$$
but individual EVPs may supply a different diagnostic
relation.
