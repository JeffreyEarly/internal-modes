---
layout: default
title: IMSurfaceGeostrophicModes
has_children: false
has_toc: false
mathjax: true
parent: Core
grand_parent: Class documentation V2
nav_order: 7
---

#  IMSurfaceGeostrophicModes

Describe surface-geostrophic boundary modes at fixed wavenumber.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMSurfaceGeostrophicModes</code></pre></div></div>

## Overview

`IMSurfaceGeostrophicModes` stores the physical boundary-value
problem for SQG boundary modes. For each horizontal wavenumber
$$k$$, the streamfunction $$\psi(z)$$ satisfies

$$
\frac{\partial}{\partial z}
\left(
\frac{f_0^2}{N^2(z)}
\frac{\partial \psi}{\partial z}
\right)
-k^2\psi=0.
$$

The active boundary is either the surface or the bottom. A surface
mode uses

$$
f_0\frac{\partial\psi}{\partial z}(z_s)=1,\qquad
f_0\frac{\partial\psi}{\partial z}(z_b)=0,
$$

while a bottom mode swaps the two endpoint conditions.

```matlab
problem = IMSurfaceGeostrophicModes.surfaceModesAtWavenumber(N2=N2,zDomain=[-4000 0],f0=1e-4,k=1e-4);
solver = IMSolverSpectral(nEVP=128);
basisSet = solver.solveSurfaceGeostrophicModes(problem);
psi = basisSet.psi(z);
```




## Topics
+ Create surface-geostrophic problems
  + [`IMSurfaceGeostrophicModes`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodes/imsurfacegeostrophicmodes.html) Create a surface-geostrophic boundary-mode problem.
  + [`bottomModesAtWavenumber`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodes/bottommodesatwavenumber.html) Create bottom SQG modes at fixed wavenumber.
  + [`surfaceModesAtWavenumber`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodes/surfacemodesatwavenumber.html) Create surface SQG modes at fixed wavenumber.
+ Summarize surface-geostrophic problems
  + [`summarize`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodes/summarize.html) Print a readable problem summary.
+ Inspect surface-geostrophic problems
  + [`N2`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodes/n2.html) Buoyancy frequency squared function.
  + [`boundary`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodes/boundary.html) Active SQG boundary, `"surface"` or `"bottom"`.
  + [`f0`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodes/f0.html) Coriolis parameter.
  + [`k`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodes/k.html) Horizontal wavenumbers.
  + [`metadata`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodes/metadata.html) Additional creation metadata.
  + [`zDomain`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodes/zdomain.html) Physical vertical domain.


---