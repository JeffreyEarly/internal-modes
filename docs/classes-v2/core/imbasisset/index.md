---
layout: default
title: IMBasisSet
has_children: false
has_toc: false
mathjax: true
parent: Core
grand_parent: Class documentation V2
nav_order: 2
---

#  IMBasisSet

Store solved EVP modes with basis-owned normalization.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMBasisSet</code></pre></div></div>

## Overview

`IMBasisSet` stores the native modes returned by a solver,
or exact analytical solution metadata supplied by a subclass.
Normalization is applied lazily when modes, spectra, or Gram matrices
are requested.

```matlab
basisSet.normalization = Normalization.geostrophic;
G = basisSet.G(linspace(-1000,0,64).');
```




## Topics
+ Create basis sets
  + [`IMBasisSet`](/internal-modes/classes-v2/core/imbasisset/imbasisset.html) Create a solved basis set.
  + [`constantStratification`](/internal-modes/classes-v2/core/imbasisset/constantstratification.html) Create an analytical constant-stratification basis set.
  + [`exponentialStratification`](/internal-modes/classes-v2/core/imbasisset/exponentialstratification.html) Create an analytical exponential-stratification basis set.
  + [`wkbApproximation`](/internal-modes/classes-v2/core/imbasisset/wkbapproximation.html) Report that WKB analytical bases are deferred.
+ Evaluate basis sets
  + [`F`](/internal-modes/classes-v2/core/imbasisset/f.html) Evaluate $$F_j(z)$$ modes on a physical grid.
  + [`G`](/internal-modes/classes-v2/core/imbasisset/g.html) Evaluate $$G_j(z)$$ modes on a physical grid.
  + [`N2`](/internal-modes/classes-v2/core/imbasisset/n2.html) Evaluate buoyancy frequency squared for this basis set.
  + [`dzLogN2`](/internal-modes/classes-v2/core/imbasisset/dzlogn2.html) Evaluate $$\partial_z\log N^2$$ for this basis set.
  + [`evaluate`](/internal-modes/classes-v2/core/imbasisset/evaluate.html) `F` or `G` on a physical grid.
  + [`evaluateAll`](/internal-modes/classes-v2/core/imbasisset/evaluateall.html) Evaluate `G` and `F` on a physical grid.
  + [`normalization`](/internal-modes/classes-v2/core/imbasisset/normalization.html) Active modal normalization.
  + [`normalizationFactors`](/internal-modes/classes-v2/core/imbasisset/normalizationfactors.html) Return factors for the requested normalization.
  + [`normalizedNativeModes`](/internal-modes/classes-v2/core/imbasisset/normalizednativemodes.html) Return native modes scaled by a normalization convention.
+ Inspect basis sets
  + [`eigenvalues`](/internal-modes/classes-v2/core/imbasisset/eigenvalues.html) Retained eigenvalues.
  + [`evp`](/internal-modes/classes-v2/core/imbasisset/evp.html) descriptor that was solved.
  + [`h`](/internal-modes/classes-v2/core/imbasisset/h.html) Equivalent depths.
  + [`index`](/internal-modes/classes-v2/core/imbasisset/index_.html) Observed and expected eigenvalue index counts.
  + [`metadata`](/internal-modes/classes-v2/core/imbasisset/metadata.html) Additional metadata.
  + [`modeNumber`](/internal-modes/classes-v2/core/imbasisset/modenumber.html) Physical mode numbers.
  + [`nativeModes`](/internal-modes/classes-v2/core/imbasisset/nativemodes.html) Native mode columns.
  + [`solver`](/internal-modes/classes-v2/core/imbasisset/solver.html) that created the native modes.
  + [`zDomain`](/internal-modes/classes-v2/core/imbasisset/zdomain.html) Physical vertical domain.
+ Analyze Gram matrices
  + [`crossSpectrum`](/internal-modes/classes-v2/core/imbasisset/crossspectrum.html) Compute a modal cross-spectrum from coefficients.
  + [`gramMatrix`](/internal-modes/classes-v2/core/imbasisset/grammatrix.html) Return the full-depth Gram matrix for a variable.
  + [`nativeTransform`](/internal-modes/classes-v2/core/imbasisset/nativetransform.html) Build a native transform.
  + [`observationProjection`](/internal-modes/classes-v2/core/imbasisset/observationprojection.html) Build an observation projection.
  + [`partialGramMatrix`](/internal-modes/classes-v2/core/imbasisset/partialgrammatrix.html) Return a partial-depth Gram matrix.
  + [`partialWindowModes`](/internal-modes/classes-v2/core/imbasisset/partialwindowmodes.html) Diagonalize partial-depth energy in existing coefficient coordinates.
  + [`spectrum`](/internal-modes/classes-v2/core/imbasisset/spectrum.html) Compute a modal spectrum from coefficients.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`emptyBoundaryWeights`](/internal-modes/classes-v2/core/imbasisset/emptyboundaryweights.html) Create an empty boundary-weight array.
  + [`scalarBoundaryWeights`](/internal-modes/classes-v2/core/imbasisset/scalarboundaryweights.html) Convert scalar endpoint weights to boundary-weight objects.


---