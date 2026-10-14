---
layout: default
title: IMGeostrophicModes
has_children: false
has_toc: false
mathjax: true
parent: Core
grand_parent: Class documentation V2
nav_order: 11
---

#  IMGeostrophicModes

Describe regular paired geostrophic eigenvalue problems.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMGeostrophicModes < IMInternalModes</code></pre></div></div>

## Overview

Numerical solves retain both fields in
$$F_\xi=-N^2G/g,\quad G_\xi=rF,\quad r=1/h.$$
Nonnegative stratification, including neutral intervals, is supported.
Spectral solves share one coordinate for both fields. Finite differences
use F on the supplied grid and G on a conservative staggered grid.
The inherited scalar coefficients describe the equivalent problem
only where N2 is positive; numerical assembly never evaluates them.
Independent zero-energy neutral-label states are not represented.

Construct through `IMInternalModes.geostrophicAPVModes` or
`IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes`.

```matlab
evp = IMInternalModes.geostrophicAPVModes(N2=N2,zDomain=[-1000 0],g0=0,gd=-0.01);
basis = IMSolverSpectral(nEVP=128).solveEVP(evp,nModes=20);
```




## Topics
+ Inspect internal-mode inner products
  + [`innerProduct`](/internal-modes/classes-v2/core/imgeostrophicmodes/innerproduct.html) Return the manuscript pairing, including companion endpoint fields.
+ Inspect mode selection
  + [`definitenessDiagnostics`](/internal-modes/classes-v2/core/imgeostrophicmodes/definitenessdiagnostics.html) Assess the regular positive F metric and endpoint energy.
  + [`modeSelectionDiagnostics`](/internal-modes/classes-v2/core/imgeostrophicmodes/modeselectiondiagnostics.html) Classify observable signed modes from the regular endpoint problem.
  + [`negativeEigenvalueBounds`](/internal-modes/classes-v2/core/imgeostrophicmodes/negativeeigenvaluebounds.html) Bound signed physical modes using endpoint inertia.
  + [`summarize`](/internal-modes/classes-v2/core/imgeostrophicmodes/summarize.html) the regular pair without evaluating scalar inverse N2.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`IMGeostrophicModes`](/internal-modes/classes-v2/core/imgeostrophicmodes/imgeostrophicmodes.html) Create a paired descriptor from a geostrophic factory configuration.
  + [`makeBasisSet`](/internal-modes/classes-v2/core/imgeostrophicmodes/makebasisset.html) Construct aligned F/G fields with stable equivalent depths.
  + [`selectModes`](/internal-modes/classes-v2/core/imgeostrophicmodes/selectmodes.html) Select using physical units without an absolute eigenvalue floor.


---