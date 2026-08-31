---
layout: default
title: IMInternalModes
has_children: false
has_toc: false
mathjax: true
parent: Core
grand_parent: Class documentation V2
nav_order: 5
---

#  IMInternalModes

Describe canonical EVPs with internal-mode interpretation.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMInternalModes < IMEigenvalueProblem</code></pre></div></div>

## Overview

`IMInternalModes` translates standard `F` and `G` internal-mode
problems into the canonical scalar EVP. The solved scalar `u` is
either `F` or `G`; the other variable is recovered diagnostically by
the relation handles `FfromGz` and `GfromFz` on the resulting
`IMInternalModesBasis`.
Internal-mode EVPs own the stratification profile `N2` and physical
vertical domain used by solvers and basis sets.

```matlab
N2 = @(z) (5.2e-3)^2*exp(2*z/1300);
evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-4000 0]);
solver = IMSolverSpectral(nEVP=128,coordinateKind="wkb");
basisSet = solver.solveEVP(evp,nModes=4);
basisSet.normalization = Normalization.geostrophic;
G = basisSet.G(linspace(-4000,0,200).');
```




## Topics
+ Create internal-mode EVPs
  + [`IMInternalModes`](/internal-modes/classes-v2/core/iminternalmodes/iminternalmodes.html) Create an internal-mode canonical EVP.
  + [`geostrophicAPVModes`](/internal-modes/classes-v2/core/iminternalmodes/geostrophicapvmodes.html) Create signed generalized-energy geostrophic APV modes.
  + [`hydrostaticFModes`](/internal-modes/classes-v2/core/iminternalmodes/hydrostaticfmodes.html) Create the hydrostatic `F` internal-mode EVP.
  + [`hydrostaticGModes`](/internal-modes/classes-v2/core/iminternalmodes/hydrostaticgmodes.html) Create the hydrostatic `G` internal-mode EVP.
  + [`meanDensityAnomalyModes`](/internal-modes/classes-v2/core/iminternalmodes/meandensityanomalymodes.html) Create generalized-energy mean-density-anomaly modes.
  + [`waveModesAtFrequency`](/internal-modes/classes-v2/core/iminternalmodes/wavemodesatfrequency.html) Create the fixed-frequency wave-mode EVP.
  + [`waveModesAtWavenumber`](/internal-modes/classes-v2/core/iminternalmodes/wavemodesatwavenumber.html) Create the fixed-wavenumber wave-mode EVP.
+ Summarize internal-mode EVPs
  + [`summarize`](/internal-modes/classes-v2/core/iminternalmodes/summarize.html) Print a readable internal-mode EVP summary.
+ Inspect internal-mode configuration
  + [`FfromGz`](/internal-modes/classes-v2/core/iminternalmodes/ffromgz.html) Diagnostic relation from `G` derivative to `F`.
  + [`GfromFz`](/internal-modes/classes-v2/core/iminternalmodes/gfromfz.html) Diagnostic relation from `F` derivative to `G`.
  + [`N2`](/internal-modes/classes-v2/core/iminternalmodes/n2.html) Buoyancy frequency squared function.
  + [`f0`](/internal-modes/classes-v2/core/iminternalmodes/f0.html) Coriolis parameter.
  + [`formulation`](/internal-modes/classes-v2/core/iminternalmodes/formulation.html) Solved physical variable, `"F"` or `"G"`.
  + [`g`](/internal-modes/classes-v2/core/iminternalmodes/g.html) Gravitational acceleration.
  + [`hFromEigenvalue`](/internal-modes/classes-v2/core/iminternalmodes/hfromeigenvalue.html) Equivalent-depth conversion function.
  + [`modeFamily`](/internal-modes/classes-v2/core/iminternalmodes/modefamily.html) Physical mode-family declaration.
+ Inspect internal-mode inner products
  + [`innerProduct`](/internal-modes/classes-v2/core/iminternalmodes/innerproduct.html) Return the `F` or `G` inner-product recipe.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`contextForSolver`](/internal-modes/classes-v2/core/iminternalmodes/contextforsolver.html) Return the internal-mode coefficient context.
  + [`makeBasisSet`](/internal-modes/classes-v2/core/iminternalmodes/makebasisset.html) Create an internal-mode basis set.


---