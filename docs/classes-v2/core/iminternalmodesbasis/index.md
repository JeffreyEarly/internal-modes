---
layout: default
title: IMInternalModesBasis
has_children: false
has_toc: false
mathjax: true
parent: Core
grand_parent: Class documentation V2
nav_order: 5
---

#  IMInternalModesBasis

Store solved internal-mode basis functions.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMInternalModesBasis < IMBasisSet</code></pre></div></div>

## Overview

`IMInternalModesBasis` evaluates the physical variables with explicit
`F(z)` and `G(z)` methods. If the EVP solves `G`, then
$$F_j=h_j G'_j.$$ If the EVP solves `F`, then
$$G_j=-gN^{-2}F'_j.$$
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
  + [`constantStratification`](/internal-modes/classes-v2/core/iminternalmodesbasis/constantstratification.html) Create an analytical constant-stratification basis set.
  + [`exponentialStratification`](/internal-modes/classes-v2/core/iminternalmodesbasis/exponentialstratification.html) Create an analytical exponential-stratification basis set.
+ Evaluate internal-mode bases
  + [`F`](/internal-modes/classes-v2/core/iminternalmodesbasis/f.html) Evaluate `F` modes.
  + [`G`](/internal-modes/classes-v2/core/iminternalmodesbasis/g.html) Evaluate `G` modes.
+ Inspect internal-mode bases
  + [`N2`](/internal-modes/classes-v2/core/iminternalmodesbasis/n2.html) Buoyancy frequency squared profile.
  + [`h`](/internal-modes/classes-v2/core/iminternalmodesbasis/h.html) Equivalent depths for the retained internal modes.


---