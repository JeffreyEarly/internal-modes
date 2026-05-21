---
layout: default
title: IMSolverWKBSpectral
has_children: false
has_toc: false
mathjax: true
parent: Solvers
grand_parent: Class documentation V2
nav_order: 3
---

#  IMSolverWKBSpectral

Solve physical-coordinate EVPs in the WKB stretched coordinate.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMSolverWKBSpectral < IMSolverSpectral</code></pre></div></div>

## Overview

The native coordinate satisfies $$dx/dz=N(z)$$. The inherited solver
applies this coordinate pullback automatically.

```matlab
solver = IMSolverWKBSpectral(N2=@(z) 1e-5*exp(z/1000), zDomain=[-1000 0], nEVP=64);
```




## Topics
+ Create solvers
  + [`IMSolverWKBSpectral`](/internal-modes/classes-v2/solvers/imsolverwkbspectral/imsolverwkbspectral.html) Create a WKB-coordinate spectral solver.


---