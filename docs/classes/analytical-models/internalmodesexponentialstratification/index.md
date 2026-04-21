---
layout: default
title: InternalModesExponentialStratification
has_children: false
has_toc: false
mathjax: true
parent: Analytical and asymptotic models
grand_parent: Class documentation
nav_order: 2
---

#  InternalModesExponentialStratification

Solve the vertical mode problem for exponential stratification.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef InternalModesExponentialStratification < InternalModesBase</code></pre></div></div>

## Overview

`InternalModesExponentialStratification` provides the closed-form
benchmark for the exponential profile used throughout the manuscript:

$$
N^2(z) = N_0^2 e^{2 z / b}
$$

together with

$$
\bar{\rho}(z) = \rho_0 \left(1 + \frac{b N_0^2}{2 g} \left(1 - e^{2 z / b}\right)\right).
$$

The resulting vertical modes are expressed in terms of Bessel
functions, matching the analytical benchmark used in the validation
sections of Early, Lelong, and Smith (2020).

```matlab
im = InternalModesExponentialStratification(N0=5.2e-3, b=1300, zIn=[-5000 0], zOut=zOut, latitude=33);
[F, G, h, omega] = im.ModesAtWavenumber(2*pi/1000);
```




## Topics
+ Create and initialize modes
  + [`InternalModesExponentialStratification`](/internal-modes/classes/analytical-models/internalmodesexponentialstratification/internalmodesexponentialstratification.html) Initialize the exponential-stratification analytical solver.
+ Inspect grids and stratification
  + [`N0`](/internal-modes/classes/analytical-models/internalmodesexponentialstratification/n0.html) Surface buoyancy frequency $$N_0$$ in radians per second.
  + [`b`](/internal-modes/classes/analytical-models/internalmodesexponentialstratification/b.html) Exponential e-folding depth $$b$$ in meters.
+ Compute modes
  + [`BottomModesAtWavenumber`](/internal-modes/classes/analytical-models/internalmodesexponentialstratification/bottommodesatwavenumber.html) Return the analytical bottom SQG mode for exponential stratification.
  + [`SurfaceModesAtWavenumber`](/internal-modes/classes/analytical-models/internalmodesexponentialstratification/surfacemodesatwavenumber.html) Return the analytical surface SQG mode for exponential stratification.
+ Inspect analytical solutions
  + [`BarotropicEquivalentDepthAtFrequency`](/internal-modes/classes/analytical-models/internalmodesexponentialstratification/barotropicequivalentdepthatfrequency.html) Estimate the barotropic equivalent depth for fixed $$\omega$$.
  + [`BarotropicEquivalentDepthAtWavenumber`](/internal-modes/classes/analytical-models/internalmodesexponentialstratification/barotropicequivalentdepthatwavenumber.html) Estimate the barotropic equivalent depth for fixed $$K$$.
  + [`FSolution`](/internal-modes/classes/analytical-models/internalmodesexponentialstratification/fsolution.html) Exact analytical $$F(z,\omega,c)$$ solution handle.
  + [`GSolution`](/internal-modes/classes/analytical-models/internalmodesexponentialstratification/gsolution.html) Exact analytical $$G(z,\omega,c)$$ solution handle.
  + [`IsStratificationExponential`](/internal-modes/classes/analytical-models/internalmodesexponentialstratification/isstratificationexponential.html) Test whether a supplied profile is close to the exponential benchmark.
  + [`N2Function`](/internal-modes/classes/analytical-models/internalmodesexponentialstratification/n2function.html) Background buoyancy-frequency function handle.
  + [`NormalizedModesForOmegaAndC`](/internal-modes/classes/analytical-models/internalmodesexponentialstratification/normalizedmodesforomegaandc.html) Evaluate and normalize analytical mode functions at `zOut`.
  + [`rhoFunction`](/internal-modes/classes/analytical-models/internalmodesexponentialstratification/rhofunction.html) Background density function handle.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`FSolutionApprox`](/internal-modes/classes/analytical-models/internalmodesexponentialstratification/fsolutionapprox.html) Approximate analytical $$F(z,\omega,c)$$ solution handle.
  + [`FindRootsInRange`](/internal-modes/classes/analytical-models/internalmodesexponentialstratification/findrootsinrange.html) Find analytical eigenvalue roots over a bounded search interval.
  + [`GSolutionApprox`](/internal-modes/classes/analytical-models/internalmodesexponentialstratification/gsolutionapprox.html) Approximate analytical $$G(z,\omega,c)$$ solution handle.
  + [`ModeFunctionsForOmegaAndC`](/internal-modes/classes/analytical-models/internalmodesexponentialstratification/modefunctionsforomegaandc.html) Select the exact or approximate analytical mode functions.
  + [`kConstantNormalizationForOmegaAndC`](/internal-modes/classes/analytical-models/internalmodesexponentialstratification/kconstantnormalizationforomegaandc.html) Return the `kConstant` normalization for the analytical exponential solution.
  + [`nInitialSearchModes`](/internal-modes/classes/analytical-models/internalmodesexponentialstratification/ninitialsearchmodes.html) Number of trial roots used when expanding the analytical search interval.
  + [`shouldApproximate`](/internal-modes/classes/analytical-models/internalmodesexponentialstratification/shouldapproximate.html) Predicate that chooses between exact and approximate Bessel forms.


---