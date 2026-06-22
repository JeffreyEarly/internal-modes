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

Describe projected surface-geostrophic boundary modes at fixed wavenumber.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMSurfaceGeostrophicModes</code></pre></div></div>

## Overview

`IMSurfaceGeostrophicModes` stores the zero-APV boundary-mode
problem for horizontal wavenumber magnitude $$k$$. The raw endpoint
modes satisfy

$$
\frac{\partial}{\partial z}
\left(
\frac{f_0^2}{N^2(z)}
\frac{\partial F}{\partial z}
\right)
-k^2F=0.
$$

A finite `g0` includes the surface buoyancy-anomaly mode; a finite
`gd` includes the bottom buoyancy-anomaly mode. Omitted or infinite
endpoint weights mean that boundary anomaly is not represented. The
solver forms the raw endpoint modes, projects them with the
boundary-energy matrix, and returns an
`IMSurfaceGeostrophicModesBasis` with `F(z)`, `G(z)`, and `h`.

The default surface anomaly includes the free-surface stretching term:

$$
\eta_0[F]=
-\frac{f_0}{N^2(0)}
\frac{\partial F}{\partial z}(0)
-\frac{f_0}{g}F(0).
$$

Use `surfaceAnomaly="noFreeSurface"` to omit the second term:

$$
\eta_0[F]=
-\frac{f_0}{N^2(0)}
\frac{\partial F}{\partial z}(0).
$$

The bottom anomaly is

$$
\eta_d[F]=
-\frac{f_0}{N^2(z_b)}
\frac{\partial F}{\partial z}(z_b).
$$

```matlab
problem = IMSurfaceGeostrophicModes.atWavenumber(N2=N2,zDomain=[-4000 0],f0=1e-4,k=1e-4,g0=-0.035);
solver = IMSolverSpectral(nEVP=128);
basisSet = solver.solveSurfaceGeostrophicModes(problem);
F = basisSet.F(z);
G = basisSet.G(z);
h = basisSet.h;
```




## Topics
+ Create surface-geostrophic problems
  + [`IMSurfaceGeostrophicModes`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodes/imsurfacegeostrophicmodes.html) Create a projected surface-geostrophic boundary-mode problem.
  + [`atWavenumber`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodes/atwavenumber.html) Create projected surface-geostrophic modes at fixed wavenumber.
+ Summarize surface-geostrophic problems
  + [`summarize`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodes/summarize.html) Print a readable problem summary.
+ Inspect surface-geostrophic problems
  + [`N2`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodes/n2.html) Buoyancy frequency squared function.
  + [`f0`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodes/f0.html) Coriolis parameter.
  + [`g`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodes/g.html) Gravitational acceleration.
  + [`g0`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodes/g0.html) Surface buoyancy-anomaly weight.
  + [`gd`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodes/gd.html) Bottom buoyancy-anomaly weight.
  + [`k`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodes/k.html) Horizontal wavenumbers.
  + [`metadata`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodes/metadata.html) Additional creation metadata.
  + [`modesPerWavenumber`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodes/modesperwavenumber.html) Return the number of projected modes for each wavenumber.
  + [`surfaceAnomaly`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodes/surfaceanomaly.html) Surface anomaly convention.
  + [`zDomain`](/internal-modes/classes-v2/core/imsurfacegeostrophicmodes/zdomain.html) Physical vertical domain.


---