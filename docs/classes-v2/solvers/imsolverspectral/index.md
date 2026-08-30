---
layout: default
title: IMSolverSpectral
has_children: false
has_toc: false
mathjax: true
parent: Solvers
grand_parent: Class documentation V2
nav_order: 2
---

#  IMSolverSpectral

Solve physical-coordinate EVPs with a Chebyshev spectral discretization.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMSolverSpectral</code></pre></div></div>

## Overview

`IMSolverSpectral` owns the numerical coordinate choice, Chebyshev
resolution, derivative matrices, and physical-coordinate pullback
rules. It is configured against an EVP or geostrophic zero-APV problem
before solving.

```matlab
evp = IMInternalModes.waveModesAtWavenumber(N2=@(z) 1e-5*ones(size(z)), zDomain=[-1000 0], k=1e-4);
solver = IMSolverSpectral(nEVP=64);
basisSet = solver.solveEVP(evp);
```




## Topics
+ Create solvers
  + [`IMSolverSpectral`](/internal-modes/classes-v2/solvers/imsolverspectral/imsolverspectral.html) Create a coordinate-aware spectral solver.
+ Inspect solvers
  + [`coordinateKind`](/internal-modes/classes-v2/solvers/imsolverspectral/coordinatekind.html) Native coordinate kind.
  + [`nEVP`](/internal-modes/classes-v2/solvers/imsolverspectral/nevp.html) Number of native EVP coefficients.
  + [`xNative`](/internal-modes/classes-v2/solvers/imsolverspectral/xnative.html) Native Lobatto grid.
  + [`zDomain`](/internal-modes/classes-v2/solvers/imsolverspectral/zdomain.html) Physical vertical domain.
  + [`zNative`](/internal-modes/classes-v2/solvers/imsolverspectral/znative.html) Physical points corresponding to `xNative`.
+ Evaluate native modes
  + [`xOfZ`](/internal-modes/classes-v2/solvers/imsolverspectral/xofz.html) Map physical coordinate to native coordinate.
  + [`zOfX`](/internal-modes/classes-v2/solvers/imsolverspectral/zofx.html) Map native coordinate to physical coordinate.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`T`](/internal-modes/classes-v2/solvers/imsolverspectral/t.html) Chebyshev basis matrix on the native grid.
  + [`Tx`](/internal-modes/classes-v2/solvers/imsolverspectral/tx.html) Native first-derivative basis matrix.
  + [`Txx`](/internal-modes/classes-v2/solvers/imsolverspectral/txx.html) Native second-derivative basis matrix.
  + [`qReference`](/internal-modes/classes-v2/solvers/imsolverspectral/qreference.html) Reference coordinate derivative $$dx/dz$$.
  + [`qzReference`](/internal-modes/classes-v2/solvers/imsolverspectral/qzreference.html) Reference physical derivative of $$dx/dz$$.
  + [`xReference`](/internal-modes/classes-v2/solvers/imsolverspectral/xreference.html) Reference native coordinate grid.
  + [`zReference`](/internal-modes/classes-v2/solvers/imsolverspectral/zreference.html) Reference physical grid for coordinate interpolation.


---