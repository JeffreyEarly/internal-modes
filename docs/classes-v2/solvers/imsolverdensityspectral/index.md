---
layout: default
title: IMSolverDensitySpectral
has_children: false
has_toc: false
mathjax: true
parent: Solvers
grand_parent: Class documentation V2
nav_order: 4
---

#  IMSolverDensitySpectral

Solve physical-coordinate EVPs in the density-like coordinate.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMSolverDensitySpectral < IMSolverSpectral</code></pre></div></div>

## Overview

The native coordinate satisfies $$dx/dz=N^2(z)$$. The inherited solver
applies this coordinate pullback automatically.

```matlab
solver = IMSolverDensitySpectral(N2=@(z) 1e-5*exp(z/1000), zDomain=[-1000 0], nEVP=64);
```




## Topics
+ Create solvers
  + [`IMSolverDensitySpectral`](/internal-modes/classes-v2/solvers/imsolverdensityspectral/imsolverdensityspectral.html) Create a density-coordinate spectral solver.


---