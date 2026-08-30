---
layout: default
title: IMInternalModesBasis
has_children: false
has_toc: false
mathjax: true
parent: Core
grand_parent: Class documentation V2
nav_order: 6
---

#  IMInternalModesBasis

Store solved internal-mode basis functions.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMInternalModesBasis < IMBasisSet</code></pre></div></div>

## Overview

`IMInternalModesBasis` evaluates the physical variables with explicit
`F(z)` and `G(z)` methods. If the EVP solves `G`, then `F(z)` is
recovered by `evp.FfromGz`, whose default relation is
$$F_j(z)=h_j\frac{\partial G_j}{\partial z}(z).$$
If the EVP solves `F`, then `G(z)` is recovered by `evp.GfromFz`.
The default inverse is hydrostatic, but wave factories install
wave-specific inverse relations.
Normalization is shared across both variables: if a rule gives scale
$$s_j$$, then both diagnostic variables for mode $$j$$ are divided by
that same factor. Standard rules are installed on the basis set, and
custom rules are added after solving with `addNormalization`.

```matlab
basisSet = solver.solveEVP(evp,nModes=4);
basisSet.normalization = Normalization.geostrophic;
F = basisSet.F(z);
G = basisSet.G(z);
```




## Topics
+ Create internal-mode bases
  + [`IMInternalModesBasis`](/internal-modes/classes-v2/core/iminternalmodesbasis/iminternalmodesbasis.html) Create an internal-mode basis set.
+ Evaluate modes
  + [`F`](/internal-modes/classes-v2/core/iminternalmodesbasis/f.html) Evaluate `F` modes.
  + [`G`](/internal-modes/classes-v2/core/iminternalmodesbasis/g.html) Evaluate `G` modes.
+ Analyze modes
  + [`crossSpectrum`](/internal-modes/classes-v2/core/iminternalmodesbasis/crossspectrum.html) Compute an internal-mode modal cross-spectrum.
  + [`gramMatrix`](/internal-modes/classes-v2/core/iminternalmodesbasis/grammatrix.html) Return a Gram matrix for `F` or `G`.
  + [`partialWindowModes`](/internal-modes/classes-v2/core/iminternalmodesbasis/partialwindowmodes.html) Diagonalize a partial-depth Gram matrix for `F` or `G`.
  + [`spectrum`](/internal-modes/classes-v2/core/iminternalmodesbasis/spectrum.html) Compute an internal-mode modal spectrum.
+ Inspect basis sets
  + [`N2`](/internal-modes/classes-v2/core/iminternalmodesbasis/n2.html) Buoyancy frequency squared profile.
  + [`h`](/internal-modes/classes-v2/core/iminternalmodesbasis/h.html) Equivalent depths for the retained internal modes.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + Normalization rules
    + [`depthNormFactor`](/internal-modes/classes-v2/core/iminternalmodesbasis/depthnormfactor.html) Return the volume-only depth normalization factor.
    + [`geostrophicNormFactor`](/internal-modes/classes-v2/core/iminternalmodesbasis/geostrophicnormfactor.html) Return the hydrostatic geostrophic normalization factor.
    + [`innerProductNormFactor`](/internal-modes/classes-v2/core/iminternalmodesbasis/innerproductnormfactor.html) Return the `F` or `G` inner-product norm factor.
    + [`maxAmplitudeNormFactor`](/internal-modes/classes-v2/core/iminternalmodesbasis/maxamplitudenormfactor.html) Return the maximum amplitude of `F` or `G`.
    + [`surfacePressureNormFactor`](/internal-modes/classes-v2/core/iminternalmodesbasis/surfacepressurenormfactor.html) Return the raw surface `F` value.
  + Gram-matrix assembly
    + [`endpointGramTerms`](/internal-modes/classes-v2/core/iminternalmodesbasis/endpointgramterms.html) Prepare rank-one endpoint terms for `F` or `G` Gram matrices.
  + Diagnostic variables
    + [`orientModeSigns`](/internal-modes/classes-v2/core/iminternalmodesbasis/orientmodesigns.html) Orient modes so the surface `F` value is positive when possible.
    + [`rawVariable`](/internal-modes/classes-v2/core/iminternalmodesbasis/rawvariable.html) Evaluate raw physical `F` or `G` modes.


---