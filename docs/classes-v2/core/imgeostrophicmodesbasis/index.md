---
layout: default
title: IMGeostrophicModesBasis
has_children: false
has_toc: false
mathjax: true
parent: Core
grand_parent: Class documentation V2
nav_order: 12
---

#  IMGeostrophicModesBasis

Store aligned geostrophic pressure and displacement modes.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMGeostrophicModesBasis < IMInternalModesBasis</code></pre></div></div>

## Overview

The regular numerical EVP solves both fields,
$$F_\xi=-N^2G/g,\qquad G_\xi=rF,\qquad r=1/h.$$
`nativeModes` remains the native `F` coefficient array. The companion
`G` coefficients are stored separately, so evaluation remains regular
wherever `N2` vanishes. Both fields share normalization and orientation.
The inverse depth is retained directly from the internal solve rather
than recovered by subtracting public physical eigenvalues.

```matlab
evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=[-1000 0],g0=0,gd=-0.01);
basis = IMSolverSpectral(nEVP=128).solveEVP(evp,nModes=10);
G = basis.G(z);
```




## Topics
+ Create internal-mode bases
  + [`IMGeostrophicModesBasis`](/internal-modes/classes-v2/core/imgeostrophicmodesbasis/imgeostrophicmodesbasis.html) Create a basis from aligned native `F` and `G` columns.
+ Evaluate modes
  + [`Gz`](/internal-modes/classes-v2/core/imgeostrophicmodesbasis/gz.html) Evaluate derivatives of the stored companion `G` modes.
+ Inspect basis sets
  + [`inverseDepth`](/internal-modes/classes-v2/core/imgeostrophicmodesbasis/inversedepth.html) Solved inverse equivalent depths, in inverse meters.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + Diagnostic variables
    + [`orientModeSigns`](/internal-modes/classes-v2/core/imgeostrophicmodesbasis/orientmodesigns.html) Apply one deterministic sign to each aligned `F`/`G` pair.
    + [`rawVariable`](/internal-modes/classes-v2/core/imgeostrophicmodesbasis/rawvariable.html) Evaluate stored unnormalized `F` or `G` modes.


---