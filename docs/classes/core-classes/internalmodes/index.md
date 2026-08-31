---
layout: default
title: InternalModes
has_children: false
has_toc: false
mathjax: true
parent: Core classes
grand_parent: Class documentation
nav_order: 1
---

#  InternalModes

Create vertical-mode solvers from gridded or analytical stratification.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef InternalModes < handle</code></pre></div></div>

## Overview

`InternalModes` is the user-facing wrapper around the concrete
solver classes. It initializes one of the numerical, analytical, or
asymptotic implementations, then forwards the common public
interface for requesting $$F_j(z)$$, $$G_j(z)$$, $$h_j$$, SQG modes, and
background-stratification diagnostics.

The wrapper follows the vertical eigenvalue problems described in
Section 2.3 of Early, Lelong, and Smith (2020):

$$
\partial_{zz} G_j = -\frac{N^2 - \omega^2}{g h_j} G_j
$$

for fixed $$\omega$$, and

$$
\partial_{zz} G_j - K^2 G_j = -\frac{N^2 - f_0^2}{g h_j} G_j
$$

for fixed $$K$$.

The primary construction path is either:

- gridded density: `InternalModes(rho, zIn, zOut, latitude, ...)`
- analytical density or `N2`: `InternalModes(rhoOrN2, [zBottom 0], zOut, latitude, ...)`

The wrapper also supports built-in benchmark profiles through the
shorthand

```matlab
im = InternalModes("constant", "wkbAdaptiveSpectral", 128);
```

and exposes quick comparison figures for those benchmark cases.

```matlab
im = InternalModes(rho, zIn, zOut, latitude, "method", "wkbAdaptiveSpectral");
[F, G, h, omega] = im.modesAtWavenumber(2*pi/1000);
psi = im.surfaceModesAtWavenumber(2*pi/1000);
```




## Topics
+ Create and initialize modes
  + [`InternalModes`](/internal-modes/classes/core-classes/internalmodes/internalmodes.html) Initialize the wrapper from a profile or a built-in benchmark case.
  + [`method`](/internal-modes/classes/core-classes/internalmodes/method.html) Name of the selected concrete solver implementation.
+ Inspect grids and stratification
  + [`Lz`](/internal-modes/classes/core-classes/internalmodes/lz.html) Total water-column depth.
  + [`N2`](/internal-modes/classes/core-classes/internalmodes/n2.html) Buoyancy frequency squared sampled on `zOut`.
  + [`N2Function`](/internal-modes/classes/core-classes/internalmodes/n2function.html) Function handle for the benchmark buoyancy-frequency profile $$N^2(z)$$.
  + [`f0`](/internal-modes/classes/core-classes/internalmodes/f0.html) Coriolis parameter at the selected latitude.
  + [`latitude`](/internal-modes/classes/core-classes/internalmodes/latitude.html) in degrees used to compute `f0`.
  + [`rho`](/internal-modes/classes/core-classes/internalmodes/rho.html) Background density sampled on `zOut`.
  + [`rho0`](/internal-modes/classes/core-classes/internalmodes/rho0.html) Reference surface density.
  + [`rhoFunction`](/internal-modes/classes/core-classes/internalmodes/rhofunction.html) Function handle for the benchmark density profile $$\bar{\rho}(z)$$.
  + [`rho_z`](/internal-modes/classes/core-classes/internalmodes/rho_z.html) First depth derivative of the background density sampled on `zOut`.
  + [`rho_zz`](/internal-modes/classes/core-classes/internalmodes/rho_zz.html) Second depth derivative of the background density sampled on `zOut`.
  + [`z`](/internal-modes/classes/core-classes/internalmodes/z.html) Output depth grid on which public profiles and modes are returned.
+ Configure normalization and boundaries
  + [`lowerBoundary`](/internal-modes/classes/core-classes/internalmodes/lowerboundary.html) Lower boundary condition at the ocean bottom.
  + [`normalization`](/internal-modes/classes/core-classes/internalmodes/normalization.html) Mode normalization convention.
  + [`upperBoundary`](/internal-modes/classes/core-classes/internalmodes/upperboundary.html) Upper boundary condition at the ocean surface.
+ Compute modes
  + [`bottomModesAtWavenumber`](/internal-modes/classes/core-classes/internalmodes/bottommodesatwavenumber.html) Return the bottom SQG mode at fixed horizontal wavenumber.
  + [`modesAtFrequency`](/internal-modes/classes/core-classes/internalmodes/modesatfrequency.html) Return vertical modes for a fixed frequency.
  + [`modesAtWavenumber`](/internal-modes/classes/core-classes/internalmodes/modesatwavenumber.html) Return vertical modes for a fixed horizontal wavenumber.
  + [`nModes`](/internal-modes/classes/core-classes/internalmodes/nmodes.html) Optional cap on the number of modes returned.
  + [`projectOntoGModesAtWavenumber`](/internal-modes/classes/core-classes/internalmodes/projectontogmodesatwavenumber.html) Project a profile onto the `G` modes at fixed horizontal wavenumber.
  + [`showLowestModesAtFrequency`](/internal-modes/classes/core-classes/internalmodes/showlowestmodesatfrequency.html) Plot the lowest resolved modes at a fixed frequency.
  + [`showLowestModesAtWavenumber`](/internal-modes/classes/core-classes/internalmodes/showlowestmodesatwavenumber.html) Plot the lowest resolved modes at a fixed horizontal wavenumber.
  + [`surfaceModesAtWavenumber`](/internal-modes/classes/core-classes/internalmodes/surfacemodesatwavenumber.html) Return the surface SQG mode at fixed horizontal wavenumber.
+ Inspect analytical solutions
  + [`ConditionNumberAsFunctionOfModeNumber`](/internal-modes/classes/core-classes/internalmodes/conditionnumberasfunctionofmodenumber.html) Compute condition number as a function of retained mode count.
  + [`NumberOfWellConditionedModes`](/internal-modes/classes/core-classes/internalmodes/numberofwellconditionedmodes.html) Estimate how many columns of a mode matrix remain numerically well conditioned.
  + [`RenormalizeForGoodConditioning`](/internal-modes/classes/core-classes/internalmodes/renormalizeforgoodconditioning.html) Renormalize a mode matrix by matching column norms to a common scale.
  + [`StratificationProfileWithName`](/internal-modes/classes/core-classes/internalmodes/stratificationprofilewithname.html) Return one of the built-in benchmark stratification profiles.
  + [`showRelativeErrorAtFrequency`](/internal-modes/classes/core-classes/internalmodes/showrelativeerroratfrequency.html) Plot benchmark relative errors for a built-in profile at fixed $$\omega$$.
  + [`showRelativeErrorAtWavenumber`](/internal-modes/classes/core-classes/internalmodes/showrelativeerroratwavenumber.html) Plot benchmark relative errors for a built-in profile at fixed $$K$$.
  + [`stratification`](/internal-modes/classes/core-classes/internalmodes/stratification.html) Name of the active built-in benchmark stratification profile.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`ConditionNumberAsFunctionOfModeNumberForModeIndices`](/internal-modes/classes/core-classes/internalmodes/conditionnumberasfunctionofmodenumberformodeindices.html) Compute mode-matrix condition numbers for selected truncation indices.
  + [`ProfileWithDoubleExponentialPlusPycnocline`](/internal-modes/classes/core-classes/internalmodes/profilewithdoubleexponentialpluspycnocline.html) Build a two-exponential pycnocline profile used by the built-in benchmarks.
  + [`ProfileWithDoubleGaussianExponentialPlusPycnocline`](/internal-modes/classes/core-classes/internalmodes/profilewithdoublegaussianexponentialpluspycnocline.html) Build a mixed Gaussian-exponential pycnocline profile used by the built-in benchmarks.
  + [`ProfileWithDoubleGaussianPlusPycnocline`](/internal-modes/classes/core-classes/internalmodes/profilewithdoublegaussianpluspycnocline.html) Build a double-Gaussian pycnocline profile used by the built-in benchmarks.
  + [`internalModes`](/internal-modes/classes/core-classes/internalmodes/internalmodes_.html) Concrete solver instance created by the wrapper.
  + [`isRunningTestCase`](/internal-modes/classes/core-classes/internalmodes/isrunningtestcase.html) Flag indicating whether the wrapper was initialized from a built-in benchmark case.
  + [`shouldShowDiagnostics`](/internal-modes/classes/core-classes/internalmodes/shouldshowdiagnostics.html) Toggle diagnostic messages from the concrete solver.


---