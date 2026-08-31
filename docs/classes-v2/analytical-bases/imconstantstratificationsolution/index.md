---
layout: default
title: IMConstantStratificationSolution
has_children: false
has_toc: false
mathjax: true
parent: Analytical bases
grand_parent: Class documentation V2
nav_order: 2
---

#  IMConstantStratificationSolution

Analytical solution family for constant stratification.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMConstantStratificationSolution < IMAnalyticalSolution</code></pre></div></div>

## Overview

`IMConstantStratificationSolution` owns the closed-form formulas for
$$N^2(z)=N_0^2$$. It can create exact internal-mode bases for the
supported canonical internal-mode EVPs and exact geostrophic zero-APV
boundary-response modes.

The public `IMInternalModes.geostrophicAPVModes` descriptor selects
an exact generalized-energy APV catalog. Positive eigendepths use
trigonometric modes, a zero eigenvalue uses the affine solution, and
negative eigendepths use hyperbolic modes. Endpoint inertia determines
whether zero, one, or two negative modes precede the optional exact
zero mode and positive branch. All finite, zero, and positive-infinite
`g0` and `gd` limits are supported under both surface conventions.
APV bases use the volume-only `depth` normalization by default and
expose roots, residuals, branch labels, endpoint inertia, `g0`, `gd`,
and `surfaceBoundary` in their metadata.

Exact geostrophic zero-APV modes use scaled hyperbolic functions with
$$m=kN_0/|f_0|$$. A two-column response solve produces unit surface
and bottom coordinates under either the free-surface response `G-F`
or the rigid-lid response `G`. The older one-boundary profile formula
used unit `f0*Fz` at its active boundary; the public exact basis instead
uses the canonical response normalization shared with the numerical
zero-APV basis.

```matlab
solution = IMConstantStratificationSolution(N0=5.2e-3,zDomain=[-5000 0]);
evp = IMInternalModes.hydrostaticGModes(N2=@(z) solution.N2(z), zDomain=solution.zDomain);
basisSet = solution.internalModes(evp,nModes=4);
exactModes = solution.geostrophicZeroAPVModesAtWavenumber(1e-4);
```




## Topics
+ Create analytical solutions
  + [`IMConstantStratificationSolution`](/internal-modes/classes-v2/analytical-bases/imconstantstratificationsolution/imconstantstratificationsolution.html) Create a constant-stratification analytical solution family.
+ Compute internal modes
  + [`internalModes`](/internal-modes/classes-v2/analytical-bases/imconstantstratificationsolution/internalmodes.html) Create an exact internal-mode basis.
+ Compute geostrophic zero-APV modes
  + [`geostrophicZeroAPVModesAtWavenumber`](/internal-modes/classes-v2/analytical-bases/imconstantstratificationsolution/geostrophiczeroapvmodesatwavenumber.html) Create exact canonical geostrophic zero-APV modes.
+ Inspect analytical solutions
  + [`N0`](/internal-modes/classes-v2/analytical-bases/imconstantstratificationsolution/n0.html) Constant buoyancy frequency $$N_0$$ in radians per second.
  + [`N2`](/internal-modes/classes-v2/analytical-bases/imconstantstratificationsolution/n2.html) Evaluate $$N^2(z)=N_0^2$$.
  + [`summarize`](/internal-modes/classes-v2/analytical-bases/imconstantstratificationsolution/summarize.html) Print a readable solution-family summary.


---