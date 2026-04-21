---
layout: default
title: InternalModesConstantStratification
has_children: false
has_toc: false
mathjax: true
parent: Analytical and asymptotic models
grand_parent: Class documentation
nav_order: 1
---

#  InternalModesConstantStratification

Solve the vertical mode problem for constant buoyancy frequency.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef InternalModesConstantStratification < InternalModesBase</code></pre></div></div>

## Overview

`InternalModesConstantStratification` provides the closed-form
solutions used throughout the package for validation and smoke
testing. In this model,

$$
N^2(z) = N_0^2
$$

and the background density is the linear profile

$$
\bar{\rho}(z) = \rho_0 - \frac{N_0^2 \rho_0}{g} z.
$$

The fixed-$$K$$ and fixed-$$\omega$$ solutions reduce to trigonometric
or hyperbolic functions with analytically known equivalent depths,
making this class the main reference implementation for checking the
numerical solvers.

```matlab
im = InternalModesConstantStratification(N0=5.2e-3, zIn=[-5000 0], zOut=zOut, latitude=33);
[F, G, h, omega] = im.ModesAtWavenumber(2*pi/1000);
```




## Topics
+ Create and initialize modes
  + [`InternalModesConstantStratification`](/internal-modes/classes/analytical-models/internalmodesconstantstratification/internalmodesconstantstratification.html) Initialize the constant-stratification analytical solver.
+ Inspect grids and stratification
  + [`N0`](/internal-modes/classes/analytical-models/internalmodesconstantstratification/n0.html) Constant buoyancy frequency $$N_0$$ in radians per second.
+ Compute modes
  + [`BottomModesAtWavenumber`](/internal-modes/classes/analytical-models/internalmodesconstantstratification/bottommodesatwavenumber.html) Return the analytical bottom SQG mode for constant stratification.
  + [`SurfaceModesAtWavenumber`](/internal-modes/classes/analytical-models/internalmodesconstantstratification/surfacemodesatwavenumber.html) Return the analytical surface SQG mode for constant stratification.
+ Inspect analytical solutions
  + [`BaroclinicModesWithEigenvalue`](/internal-modes/classes/analytical-models/internalmodesconstantstratification/baroclinicmodeswitheigenvalue.html) Evaluate the analytical baroclinic mode shapes for given eigenvalues.
  + [`BarotropicModeAtFrequency`](/internal-modes/classes/analytical-models/internalmodesconstantstratification/barotropicmodeatfrequency.html) Return the analytical barotropic mode branch for fixed $$\omega$$.
  + [`BarotropicModeAtWavenumber`](/internal-modes/classes/analytical-models/internalmodesconstantstratification/barotropicmodeatwavenumber.html) Return the analytical barotropic mode branch for fixed $$K$$.
  + [`BuoyancyFrequencyFromConstantStratification`](/internal-modes/classes/analytical-models/internalmodesconstantstratification/buoyancyfrequencyfromconstantstratification.html) Estimate `N0` and `rho0` from a constant-stratification profile.
  + [`IsStratificationConstant`](/internal-modes/classes/analytical-models/internalmodesconstantstratification/isstratificationconstant.html) Test whether a supplied profile is close to constant stratification.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`BaroclinicModeNormalization`](/internal-modes/classes/analytical-models/internalmodesconstantstratification/baroclinicmodenormalization.html) Return the requested analytical normalization factor for baroclinic modes.
  + [`BarotropicMode`](/internal-modes/classes/analytical-models/internalmodesconstantstratification/barotropicmode.html) Evaluate a chosen analytical barotropic branch.
  + [`BarotropicModeNormalization`](/internal-modes/classes/analytical-models/internalmodesconstantstratification/barotropicmodenormalization.html) Return the requested analytical normalization factor for the barotropic branch.


---