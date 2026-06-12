---
layout: default
title: selectModes
parent: IMIndexPolicy
grand_parent: Core
nav_order: 13
mathjax: true
---

#  selectModes

Select and label retained modes according to the index policy.


---

## Declaration
```matlab
 selection = selectModes(policy,eigenvalues,nModes,context)
```
## Parameters
+ `eigenvalues`  candidate eigenvalues
+ `nModes`  number of modes to retain
+ `context`  solver or EVP context

## Returns
+ `selection`  structure with `sortIndex`, `modeNumber`, and `index`

## Discussion

  The returned `modeNumber` uses `-1` for a surface boundary
  mode, `-2` for a bottom boundary mode, `0` for a true null
  mode with $$F_0(z)=1$$ and $$G_0(z)=0$$, and positive labels
  for interior baroclinic modes.
