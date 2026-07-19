---
layout: default
title: IMBasisSet
has_children: false
has_toc: false
mathjax: true
parent: Core
grand_parent: Class documentation V2
nav_order: 3
---

#  IMBasisSet

Store solved scalar canonical EVP modes.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMBasisSet</code></pre></div></div>

## Overview

`IMBasisSet` stores the scalar modes selected from a canonical EVP and
evaluates the solved variable `u` and derivative `uz`. Modal
normalization is applied lazily when values, Gram matrices, or spectra
are requested. For each retained mode,
$$u_j^{\mathrm{out}}(z)=u_j^{\mathrm{raw}}(z)/s_j,$$
where the scale factor $$s_j$$ is supplied by a basis-set
normalization rule. Custom rules are added after solving with
`addNormalization`.

```matlab
basisSet = solver.solveEVP(evp,nModes=4);
basisSet = basisSet.addNormalization("constantScaled", @(basisSet,j) C*basisSet.innerProductNormFactor(j));
basisSet.normalization = "unity";
u = basisSet.u(z);
factors = basisSet.normalizationFactors("unity");
```




## Topics
+ Create basis sets
  + [`IMBasisSet`](/internal-modes/classes-v2/core/imbasisset/imbasisset.html) Create a solved scalar basis set.
+ Evaluate modes
  + [`u`](/internal-modes/classes-v2/core/imbasisset/u.html) Evaluate solved scalar modes.
  + [`uz`](/internal-modes/classes-v2/core/imbasisset/uz.html) Evaluate solved scalar vertical derivatives.
+ Normalize modes
  + [`addNormalization`](/internal-modes/classes-v2/core/imbasisset/addnormalization.html) Add a named normalization rule.
  + [`normalization`](/internal-modes/classes-v2/core/imbasisset/normalization.html) Active normalization rule name.
  + [`normalizationFactors`](/internal-modes/classes-v2/core/imbasisset/normalizationfactors.html) Return scale factors for a normalization rule.
  + [`normalizationNames`](/internal-modes/classes-v2/core/imbasisset/normalizationnames.html) Return installed normalization rule names.
+ Analyze modes
  + [`crossSpectrum`](/internal-modes/classes-v2/core/imbasisset/crossspectrum.html) Compute a scalar modal cross-spectrum.
  + [`gramMatrix`](/internal-modes/classes-v2/core/imbasisset/grammatrix.html) Return a scalar Gram matrix.
  + [`partialWindowModes`](/internal-modes/classes-v2/core/imbasisset/partialwindowmodes.html) Diagonalize a partial scalar Gram matrix.
  + [`spectrum`](/internal-modes/classes-v2/core/imbasisset/spectrum.html) Compute a scalar modal spectrum.
+ Build discrete transforms
  + [`discreteTransform`](/internal-modes/classes-v2/core/imbasisset/discretetransform.html) Build a scalar Galerkin transform on fixed sample points.
  + [`pointsFromModeRoots`](/internal-modes/classes-v2/core/imbasisset/pointsfrommoderoots.html) Return physical endpoints and roots of the next selected mode.
  + [`quadratureWeightsForPoints`](/internal-modes/classes-v2/core/imbasisset/quadratureweightsforpoints.html) Find quadrature weights for fixed physical sample points.
+ Inspect basis sets
  + [`eigenvalues`](/internal-modes/classes-v2/core/imbasisset/eigenvalues.html) Retained eigenvalues.
  + [`evp`](/internal-modes/classes-v2/core/imbasisset/evp.html) descriptor that was solved.
  + [`metadata`](/internal-modes/classes-v2/core/imbasisset/metadata.html) Additional metadata.
  + [`modeNumber`](/internal-modes/classes-v2/core/imbasisset/modenumber.html) Retained-mode labels.
  + [`modeSelectionDiagnostics`](/internal-modes/classes-v2/core/imbasisset/modeselectiondiagnostics.html) Mode-selection diagnostics.
  + [`zDomain`](/internal-modes/classes-v2/core/imbasisset/zdomain.html) Physical vertical domain.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + Gram-matrix assembly
    + [`endpointGramTerms`](/internal-modes/classes-v2/core/imbasisset/endpointgramterms.html) Prepare rank-one endpoint terms for scalar Gram matrices.
  + Normalization rules
    + [`innerProductNormFactor`](/internal-modes/classes-v2/core/imbasisset/innerproductnormfactor.html) Return the scalar inner-product norm factor.
    + [`maxAmplitudeNormFactor`](/internal-modes/classes-v2/core/imbasisset/maxamplitudenormfactor.html) Return the maximum scalar amplitude.
  + Native representation
    + [`nativeModes`](/internal-modes/classes-v2/core/imbasisset/nativemodes.html) Native mode columns.
    + [`orientModeSigns`](/internal-modes/classes-v2/core/imbasisset/orientmodesigns.html) Orient scalar modes with a deterministic sign convention.
    + [`rawU`](/internal-modes/classes-v2/core/imbasisset/rawu.html) Evaluate raw solved scalar modes.
    + [`rawUz`](/internal-modes/classes-v2/core/imbasisset/rawuz.html) Evaluate raw solved scalar derivatives.
    + [`solver`](/internal-modes/classes-v2/core/imbasisset/solver.html) that created the native modes.


---