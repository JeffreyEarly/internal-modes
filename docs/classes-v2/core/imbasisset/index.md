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

`IMBasisSet` stores native coefficient columns returned by a solver
and evaluates the solved scalar variable `u` and its derivative. Modal
normalization is applied lazily when values or Gram matrices are
requested. For each mode,
$$u_j^{\mathrm{out}}(z)=u_j^{\mathrm{raw}}(z)/s_j,$$
where the scale factor $$s_j$$ is supplied by the active
normalization rule.

```matlab
basisSet = solver.solveEVP(evp,nModes=4);
basisSet.normalization = Normalization.unity;
u = basisSet.u(z);
factors = basisSet.normalizationFactors(Normalization.unity);
```




## Topics
+ Create basis sets
  + [`IMBasisSet`](/internal-modes/classes-v2/core/imbasisset/imbasisset.html) Create a solved scalar basis set.
  + [`constantStratification`](/internal-modes/classes-v2/core/imbasisset/constantstratification.html) Create an analytical constant-stratification basis set.
  + [`exponentialStratification`](/internal-modes/classes-v2/core/imbasisset/exponentialstratification.html) Create an analytical exponential-stratification basis set.
+ Evaluate basis sets
  + [`N2`](/internal-modes/classes-v2/core/imbasisset/n2.html) Evaluate buoyancy frequency squared.
  + [`dzLogN2`](/internal-modes/classes-v2/core/imbasisset/dzlogn2.html) Evaluate the vertical derivative of `log(N2)`.
  + [`evaluate`](/internal-modes/classes-v2/core/imbasisset/evaluate.html) the scalar variable.
  + [`evaluateAll`](/internal-modes/classes-v2/core/imbasisset/evaluateall.html) Evaluate all scalar fields.
  + [`normalization`](/internal-modes/classes-v2/core/imbasisset/normalization.html) Active modal normalization.
  + [`normalizationFactors`](/internal-modes/classes-v2/core/imbasisset/normalizationfactors.html) Return factors for a normalization convention.
  + [`normalizedNativeModes`](/internal-modes/classes-v2/core/imbasisset/normalizednativemodes.html) Return native modes scaled by a normalization.
  + [`u`](/internal-modes/classes-v2/core/imbasisset/u.html) Evaluate solved scalar modes.
  + [`uz`](/internal-modes/classes-v2/core/imbasisset/uz.html) Evaluate solved scalar vertical derivatives.
+ Analyze Gram matrices
  + [`crossSpectrum`](/internal-modes/classes-v2/core/imbasisset/crossspectrum.html) Compute a scalar modal cross-spectrum.
  + [`gramMatrix`](/internal-modes/classes-v2/core/imbasisset/grammatrix.html) Return the full-domain scalar Gram matrix.
  + [`partialGramMatrix`](/internal-modes/classes-v2/core/imbasisset/partialgrammatrix.html) Return a partial-domain scalar Gram matrix.
  + [`partialWindowModes`](/internal-modes/classes-v2/core/imbasisset/partialwindowmodes.html) Diagonalize a partial scalar Gram matrix.
  + [`spectrum`](/internal-modes/classes-v2/core/imbasisset/spectrum.html) Compute a scalar modal spectrum.
+ Inspect basis sets
  + [`eigenvalues`](/internal-modes/classes-v2/core/imbasisset/eigenvalues.html) Retained eigenvalues.
  + [`evp`](/internal-modes/classes-v2/core/imbasisset/evp.html) descriptor that was solved.
  + [`h`](/internal-modes/classes-v2/core/imbasisset/h.html) Equivalent depths.
  + [`index`](/internal-modes/classes-v2/core/imbasisset/index_.html) Diagnostic index information from mode selection.
  + [`metadata`](/internal-modes/classes-v2/core/imbasisset/metadata.html) Additional metadata.
  + [`modeNumber`](/internal-modes/classes-v2/core/imbasisset/modenumber.html) Physical mode numbers.
  + [`nativeModes`](/internal-modes/classes-v2/core/imbasisset/nativemodes.html) Native mode columns.
  + [`solver`](/internal-modes/classes-v2/core/imbasisset/solver.html) that created the native modes.
  + [`zDomain`](/internal-modes/classes-v2/core/imbasisset/zdomain.html) Physical vertical domain.


---