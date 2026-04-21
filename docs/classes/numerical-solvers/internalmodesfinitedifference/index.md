---
layout: default
title: InternalModesFiniteDifference
has_children: false
has_toc: false
mathjax: true
parent: Numerical solvers
grand_parent: Class documentation
nav_order: 5
---

#  InternalModesFiniteDifference

Solve the vertical eigenvalue problems with finite-difference matrices.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef InternalModesFiniteDifference < InternalModesBase</code></pre></div></div>

## Overview

`InternalModesFiniteDifference` discretizes the fixed-$$K$$ and
fixed-$$\omega$$ eigenvalue problems directly on the supplied depth
grid. It uses arbitrary-order finite-difference stencils generated
from the Fornberg weights algorithm and optionally interpolates the
resulting modes onto a separate output grid.

This class is mainly useful as a baseline against the spectral
methods discussed in Section 3 of Early, Lelong, and Smith (2020),
where the manuscript explains why finite differencing becomes less
attractive when many frequencies or wavenumbers are required.

```matlab
im = InternalModesFiniteDifference(rho=rho, zIn=zIn, zOut=zOut, latitude=latitude, orderOfAccuracy=4);
[F, G, h, omega] = im.modesAtWavenumber(2*pi/1000);
```




## Topics
+ Create and initialize modes
  + [`InternalModesFiniteDifference`](/internal-modes/classes/numerical-solvers/internalmodesfinitedifference/internalmodesfinitedifference.html) Initialize the finite-difference solver.
  + [`orderOfAccuracy`](/internal-modes/classes/numerical-solvers/internalmodesfinitedifference/orderofaccuracy.html) Formal order of accuracy for the finite-difference stencils.
+ Compute modes
  + [`bottomModesAtWavenumber`](/internal-modes/classes/numerical-solvers/internalmodesfinitedifference/bottommodesatwavenumber.html) Return the bottom SQG mode at fixed horizontal wavenumber.
  + [`boundaryModesAtWavenumber`](/internal-modes/classes/numerical-solvers/internalmodesfinitedifference/boundarymodesatwavenumber.html) Return either the surface or bottom SQG mode at fixed horizontal wavenumber.
  + [`surfaceModesAtWavenumber`](/internal-modes/classes/numerical-solvers/internalmodesfinitedifference/surfacemodesatwavenumber.html) Return the surface SQG mode at fixed horizontal wavenumber.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`ApplyBoundaryConditions`](/internal-modes/classes/numerical-solvers/internalmodesfinitedifference/applyboundaryconditions.html) Apply the current lower and upper boundary conditions to an EVP pair.
  + [`Diff1`](/internal-modes/classes/numerical-solvers/internalmodesfinitedifference/diff1.html) First-derivative finite-difference matrix.
  + [`Diff2`](/internal-modes/classes/numerical-solvers/internalmodesfinitedifference/diff2.html) Second-derivative finite-difference matrix.
  + [`FiniteDifferenceMatrix`](/internal-modes/classes/numerical-solvers/internalmodesfinitedifference/finitedifferencematrix.html) Build a finite-difference differentiation matrix with boundary stencils.
  + [`InitializeOutputTransformation`](/internal-modes/classes/numerical-solvers/internalmodesfinitedifference/initializeoutputtransformation.html) Prepare the interpolation from the differentiation grid to the public output grid.
  + [`ModesFromGEP`](/internal-modes/classes/numerical-solvers/internalmodesfinitedifference/modesfromgep.html) Solve a generalized EVP and map its modes onto the public grid.
  + [`N2_z_diff`](/internal-modes/classes/numerical-solvers/internalmodesfinitedifference/n2_z_diff.html) Buoyancy frequency squared sampled on `z_diff`.
  + [`NormalizeModes`](/internal-modes/classes/numerical-solvers/internalmodesfinitedifference/normalizemodes.html) Normalize finite-difference modes using the active convention.
  + [`T_out`](/internal-modes/classes/numerical-solvers/internalmodesfinitedifference/t_out.html) Transformation from `z_diff` functions to the public output grid.
  + [`n`](/internal-modes/classes/numerical-solvers/internalmodesfinitedifference/n.html) Number of differentiation-grid points.
  + [`rho_z_diff`](/internal-modes/classes/numerical-solvers/internalmodesfinitedifference/rho_z_diff.html) Background density sampled on `z_diff`.
  + [`weights`](/internal-modes/classes/numerical-solvers/internalmodesfinitedifference/weights.html) Return Fornberg finite-difference weights for one stencil location.
  + [`z_diff`](/internal-modes/classes/numerical-solvers/internalmodesfinitedifference/z_diff.html) Depth grid used for differentiation and the generalized EVP.


---