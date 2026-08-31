---
layout: default
title: GfromFz
parent: IMInternalModes
grand_parent: Core
nav_order: 2
mathjax: true
---

#  GfromFz

Diagnostic relation from `F` derivative to `G`.


---

## Discussion

  `GfromFz` has signature `G = GfromFz(z,dFdz,h,ctx)`. The default
  relation is the hydrostatic inverse
  $$G_j(z)=-\frac{g}{N^2(z)}
  \frac{\partial F_j}{\partial z}(z).$$
  Wave-mode factories install relation handles with the
  appropriate wave correction factors.
