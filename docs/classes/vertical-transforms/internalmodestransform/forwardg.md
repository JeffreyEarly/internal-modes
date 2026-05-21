---
layout: default
title: forwardG
parent: InternalModesTransform
grand_parent: Classes
nav_order: 14
mathjax: true
---

#  forwardG

Forward G projection matrix from samples to modal coefficients.


---

## Description
Real valued property with dimensions $$(modeG,zIndex)$$ and no units.

## Discussion

  For rigid-lid G modes, zero endpoint columns may be present so the
  matrix remains compatible with fields sampled on the full `z`
  grid while the projection itself uses the active interior rows.
