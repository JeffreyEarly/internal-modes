---
layout: default
title: IMMeanDensityAnomalyModesBasis
has_children: false
has_toc: false
mathjax: true
parent: Core
grand_parent: Class documentation V2
nav_order: 8
---

#  IMMeanDensityAnomalyModesBasis

Store aligned mean-density-anomaly `F` and `G` modes.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMMeanDensityAnomalyModesBasis < IMInternalModesBasis</code></pre></div></div>

## Overview

The EVP solves `G`. This basis stores the companion pressure shapes
obtained by surface-referenced integration of the native spectral or
finite-difference representation:

$$
F_j(z)=\frac{1}{g}\int_z^{z_s}N^2(z')G_j(z')\,dz'.
$$

The default `unity` normalization divides each aligned `F`/`G` pair
by the positive magnitude of its signed generalized-energy norm.
`signatures(j)` records the remaining sign, so the normalized Gram
matrix is `diag(signatures)`.

Discrete transforms project `G` and synthesize both `G` and the
diagnostic `F`. Omitted transform variables therefore select `G`.




## Topics
+ Create internal-mode bases
  + [`IMMeanDensityAnomalyModesBasis`](/internal-modes/classes-v2/core/immeandensityanomalymodesbasis/immeandensityanomalymodesbasis.html) Create an aligned mean-density-anomaly basis set.
+ Evaluate modes
  + [`Gz`](/internal-modes/classes-v2/core/immeandensityanomalymodesbasis/gz.html) Evaluate vertical derivatives of the solved `G` modes.
+ Analyze modes
  + [`signatures`](/internal-modes/classes-v2/core/immeandensityanomalymodesbasis/signatures.html) Signs of the normalized generalized-energy mode norms.
+ Inspect basis sets
  + [`activeEndpoints`](/internal-modes/classes-v2/core/immeandensityanomalymodesbasis/activeendpoints.html) Active endpoint identities in canonical surface-bottom order.
  + [`g0`](/internal-modes/classes-v2/core/immeandensityanomalymodesbasis/g0.html) Surface generalized-energy acceleration.
  + [`gd`](/internal-modes/classes-v2/core/immeandensityanomalymodesbasis/gd.html) Bottom generalized-energy acceleration.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + Diagnostic variables
    + [`orientModeSigns`](/internal-modes/classes-v2/core/immeandensityanomalymodesbasis/orientmodesigns.html) Apply one deterministic sign to each aligned `F`/`G` pair.
    + [`rawVariable`](/internal-modes/classes-v2/core/immeandensityanomalymodesbasis/rawvariable.html) Evaluate unnormalized aligned `F` or `G` modes.


---