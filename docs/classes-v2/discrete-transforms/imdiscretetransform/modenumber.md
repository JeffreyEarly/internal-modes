---
layout: default
title: modeNumber
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 10
mathjax: true
---

#  modeNumber

Physical labels for the retained modal rows and columns.


---

## Discussion

`modeNumber` is a $$1\times n_m$$ row vector. Entry $$j$$ labels
column $$j$$ of `inverseMatrix`, row $$j$$ of `forwardMatrix`, and
row $$j$$ of arrays returned by `transformForward`. These are mode
labels, not MATLAB array indices, and may include negative or zero
values when the source basis contains such modes.
