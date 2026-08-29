---
layout: default
title: IMExponentialStratificationSolution
has_children: false
has_toc: false
mathjax: true
parent: Analytical bases
grand_parent: Class documentation V2
nav_order: 3
---

#  IMExponentialStratificationSolution

Analytical solution family for exponential stratification.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMExponentialStratificationSolution < IMAnalyticalSolution</code></pre></div></div>

## Overview

`IMExponentialStratificationSolution` owns the closed-form formulas for
$$N^2(z)=N_0^2\exp(2z/b)$$ on domains with the surface at $$z=0$$. It
can create exact internal-mode bases for supported rigid-bottom
internal-mode EVPs, generalized-energy geostrophic APV EVPs, and
exact surface or bottom SQG boundary modes.

The generalized-energy APV branch recognizes a hydrostatic `F` EVP
named `"geostrophicAPVModes"` with canonical coefficients
$$p=1/N^2$$, $$q=0$$, and $$r=1/g$$. Its parameter struct must contain
the signed endpoint accelerations `g0` and `gd`, together with
`surfaceBoundary="freeSurface"` or `"rigidLid"`. Finite, zero, and
positive-infinite endpoint values are supported. The returned exact
basis is ordered by $$1/h$$ and may contain negative modes, an exact
zero-eigenvalue mode represented by `h=Inf`, and positive modes.

```matlab
solution = IMExponentialStratificationSolution(N0=5.2e-3,b=1300,zDomain=[-5000 0]);
evp = IMInternalModes.hydrostaticGModes(N2=@(z) solution.N2(z), zDomain=solution.zDomain);
basisSet = solution.internalModes(evp,nModes=4);
```




## Topics
+ Create analytical solutions
  + [`IMExponentialStratificationSolution`](/internal-modes/classes-v2/analytical-bases/imexponentialstratificationsolution/imexponentialstratificationsolution.html) Create an exponential-stratification analytical solution family.
+ Compute internal modes
  + [`internalModeAvailability`](/internal-modes/classes-v2/analytical-bases/imexponentialstratificationsolution/internalmodeavailability.html) Report whether exact internal modes are available.
  + [`internalModes`](/internal-modes/classes-v2/analytical-bases/imexponentialstratificationsolution/internalmodes.html) Create an exact internal-mode basis.
+ Compute SQG modes
  + [`sqgAvailability`](/internal-modes/classes-v2/analytical-bases/imexponentialstratificationsolution/sqgavailability.html) Report whether exact SQG modes are available.
  + [`sqgModesAtWavenumber`](/internal-modes/classes-v2/analytical-bases/imexponentialstratificationsolution/sqgmodesatwavenumber.html) Create exact SQG boundary modes at fixed wavenumber.
+ Inspect analytical solutions
  + [`N0`](/internal-modes/classes-v2/analytical-bases/imexponentialstratificationsolution/n0.html) Surface buoyancy frequency $$N_0$$ in radians per second.
  + [`N2`](/internal-modes/classes-v2/analytical-bases/imexponentialstratificationsolution/n2.html) Evaluate $$N^2(z)=N_0^2\exp(2z/b)$$.
  + [`b`](/internal-modes/classes-v2/analytical-bases/imexponentialstratificationsolution/b.html) Exponential e-folding depth $$b$$ in meters.
  + [`summarize`](/internal-modes/classes-v2/analytical-bases/imexponentialstratificationsolution/summarize.html) Print a readable solution-family summary.


---