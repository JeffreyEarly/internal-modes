---
layout: default
title: InternalModesBase
has_children: false
has_toc: false
mathjax: true
parent: Core classes
grand_parent: Class documentation
nav_order: 2
---

#  InternalModesBase

Define the shared contract for all internal-mode solvers.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef (Abstract) InternalModesBase < handle</code></pre></div></div>

## Overview

`InternalModesBase` is the abstract base class for the concrete
solver implementations. It stores the common physical parameters,
normalization choices, and boundary-condition state, and it defines
the shared public interface for requesting modes at fixed wavenumber
or fixed frequency.

Following Section 2.3 of Early, Lelong, and Smith (2020), the
package centers on the vertical structure functions `F_j(z)` and
`G_j(z)` with equivalent depths `h_j`, connected by

$$
(N^2 - \omega^2) G_j = -g \, \partial_z F_j,
\qquad
F_j = h_j \, \partial_z G_j.
$$

Concrete subclasses solve either the fixed-`K` or fixed-`\omega`
eigenvalue problems using spectral, finite-difference, analytical,
or WKB-based formulations, but they all expose the same shared state
documented here.




## Topics
+ Inspect grids and stratification
  + [`Lz`](/internal-modes/classes/core-classes/internalmodesbase/lz.html) Total water-column depth `D = zMax - zMin`.
  + [`N2`](/internal-modes/classes/core-classes/internalmodesbase/n2.html) Buoyancy frequency squared sampled on `zOut`.
  + [`f0`](/internal-modes/classes/core-classes/internalmodesbase/f0.html) Coriolis parameter at the selected latitude.
  + [`latitude`](/internal-modes/classes/core-classes/internalmodesbase/latitude.html) in degrees used to compute `f0`.
  + [`rho`](/internal-modes/classes/core-classes/internalmodesbase/rho.html) Background density sampled on `zOut`.
  + [`rho0`](/internal-modes/classes/core-classes/internalmodesbase/rho0.html) Reference surface density.
  + [`rho_z`](/internal-modes/classes/core-classes/internalmodesbase/rho_z.html) First depth derivative of the background density sampled on `zOut`.
  + [`rho_zz`](/internal-modes/classes/core-classes/internalmodesbase/rho_zz.html) Second depth derivative of the background density sampled on `zOut`.
  + [`rotationRate`](/internal-modes/classes/core-classes/internalmodesbase/rotationrate.html) Planetary rotation rate in radians per second.
  + [`z`](/internal-modes/classes/core-classes/internalmodesbase/z.html) Output depth grid on which public profiles and modes are returned.
  + [`zDomain`](/internal-modes/classes/core-classes/internalmodesbase/zdomain.html) Two-element depth domain `[zMin zMax]`.
  + [`zMax`](/internal-modes/classes/core-classes/internalmodesbase/zmax.html) Maximum depth in `zDomain`.
  + [`zMin`](/internal-modes/classes/core-classes/internalmodesbase/zmin.html) Minimum depth in `zDomain`.
+ Configure normalization and boundaries
  + [`lowerBoundary`](/internal-modes/classes/core-classes/internalmodesbase/lowerboundary.html) Lower boundary condition at the ocean bottom.
  + [`normalization`](/internal-modes/classes/core-classes/internalmodesbase/normalization.html) Mode normalization convention.
  + [`upperBoundary`](/internal-modes/classes/core-classes/internalmodesbase/upperboundary.html) Upper boundary condition at the ocean surface.
+ Compute modes
  + [`ModesAtFrequency`](/internal-modes/classes/core-classes/internalmodesbase/modesatfrequency.html) Return vertical modes for a fixed frequency.
  + [`ModesAtWavenumber`](/internal-modes/classes/core-classes/internalmodesbase/modesatwavenumber.html) Return vertical modes for a fixed horizontal wavenumber.
  + [`nModes`](/internal-modes/classes/core-classes/internalmodesbase/nmodes.html) Optional cap on the number of modes returned.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`gridFrequency`](/internal-modes/classes/core-classes/internalmodesbase/gridfrequency.html) Last fixed frequency used to build an adaptive grid, when applicable.
  + [`requiresMonotonicDensity`](/internal-modes/classes/core-classes/internalmodesbase/requiresmonotonicdensity.html) Flag indicating whether the concrete solver requires monotonic density.
  + [`shouldShowDiagnostics`](/internal-modes/classes/core-classes/internalmodesbase/shouldshowdiagnostics.html) Toggle diagnostic messages from the active solver.


---