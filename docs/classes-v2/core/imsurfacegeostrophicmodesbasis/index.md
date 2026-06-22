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

Store projected surface-geostrophic boundary modes.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMSurfaceGeostrophicModesBasis</code></pre></div></div>

## Overview

`IMSurfaceGeostrophicModesBasis` evaluates the projected zero-APV
boundary modes computed by a numerical solver from an
`IMSurfaceGeostrophicModes` problem. The solved mode amplitude is
`F(z)`. The diagnostic displacement variable is

$$
G(z)=
-\frac{g}{N^2(z)}
\frac{\partial F}{\partial z}(z).
$$

The equivalent boundary depths satisfy

$$
h_0^a=2k^2\gamma_a,
$$

where $$\gamma_a$$ are the eigenvalues of the boundary-energy
matrix used to project the raw endpoint modes.

```matlab
basisSet = solver.solveSurfaceGeostrophicModes(problem);
F = basisSet.F(z);
G = basisSet.G(z);
h = basisSet.h;
```




## Topics
+ Evaluate surface-geostrophic modes
  + [`F`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/f.html) Evaluate projected SQG streamfunction modes.
  + [`G`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/g.html) Evaluate diagnostic SQG displacement modes.
  + [`IMSurfaceGeostrophicModesBasis`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/imsurfacegeostrophicmodesbasis.html) Create a solved surface-geostrophic basis.
+ Inspect surface-geostrophic modes
  + [`N2`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/n2.html) Buoyancy frequency squared function.
  + [`energyEigenvalues`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/energyeigenvalues.html) Boundary-energy eigenvalues.
  + [`g`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/g_.html) Gravitational acceleration.
  + [`g0`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/g0.html) Surface buoyancy-anomaly weight.
  + [`gd`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/gd.html) Bottom buoyancy-anomaly weight.
  + [`h`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/h.html) Equivalent boundary depths.
  + [`k`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/k.html) Mode-aligned horizontal wavenumbers.
  + [`metadata`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/metadata.html) Additional creation metadata.
  + [`mixingCoefficients`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/mixingcoefficients.html) Raw endpoint-mode mixing coefficients.
  + [`modeNumber`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/modenumber.html) Projected mode labels within each wavenumber.
  + [`nativeModes`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/nativemodes.html) Native projected mode columns before interpolation.
  + [`problem`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/problem.html) Surface-geostrophic problem descriptor.
  + [`solver`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/solver.html) Numerical solver used to compute the native mode columns.
  + [`summarize`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/summarize.html) Print a readable surface-geostrophic basis summary.
  + [`surfaceAnomaly`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/surfaceanomaly.html) Surface anomaly convention.
  + [`zDomain`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodesbasis/zdomain.html) Physical vertical domain.


---