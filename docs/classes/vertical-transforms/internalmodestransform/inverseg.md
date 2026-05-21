---
layout: default
title: inverseG
parent: InternalModesTransform
grand_parent: Classes
nav_order: 24
mathjax: true
---

#  inverseG

Inverse G reconstruction matrix from coefficients to samples.


---

## Description
Real valued property with dimensions $$(zIndex,modeG)$$ and no units.

## Discussion

  For rigid-lid G modes, zero endpoint rows may be present so the
  matrix reconstructs onto the full `z` grid.
