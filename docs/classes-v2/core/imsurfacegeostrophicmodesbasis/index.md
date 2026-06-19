---
layout: default
title: IMSurfaceGeostrophicModesBasis
has_children: false
has_toc: false
mathjax: true
parent: Core
grand_parent: Class documentation V2
nav_order: 8
---

#  IMSurfaceGeostrophicModesBasis

Store solved surface-geostrophic boundary modes.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMSurfaceGeostrophicModesBasis</code></pre></div></div>

## Overview

`IMSurfaceGeostrophicModesBasis` evaluates the SQG streamfunction
modes computed by a numerical solver from an
`IMSurfaceGeostrophicModes` problem.

```matlab
basisSet = solver.solveSurfaceGeostrophicModes(problem);
psi = basisSet.psi(z);
psiz = basisSet.psiz(z);
```




## Topics
+ Evaluate surface-geostrophic modes
  + [`IMSurfaceGeostrophicModesBasis`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/imsurfacegeostrophicmodesbasis.html) Create a solved surface-geostrophic basis.
  + [`psi`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/psi.html) Evaluate SQG streamfunction modes.
  + [`psiz`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/psiz.html) Evaluate vertical derivatives of SQG streamfunction modes.
+ Inspect surface-geostrophic modes
  + [`N2`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/n2.html) Buoyancy frequency squared function.
  + [`boundary`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/boundary.html) Active SQG boundary, `"surface"` or `"bottom"`.
  + [`k`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/k.html) Horizontal wavenumbers.
  + [`metadata`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/metadata.html) Additional creation metadata.
  + [`nativeModes`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/nativemodes.html) Native mode columns before interpolation.
  + [`problem`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/problem.html) Surface-geostrophic problem descriptor.
  + [`solver`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/solver.html) Numerical solver used to compute the native mode columns.
  + [`summarize`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/summarize.html) Print a readable surface-geostrophic basis summary.
  + [`zDomain`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/zdomain.html) Physical vertical domain.


---