---
layout: default
title: InternalModesWKB
has_children: false
has_toc: false
mathjax: true
parent: Analytical and asymptotic models
grand_parent: Class documentation
nav_order: 3
---

#  InternalModesWKB

Compute WKB mode approximations from a spectrally resolved stratification.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef InternalModesWKB < InternalModesSpectral</code></pre></div></div>

## Overview

`InternalModesWKB` uses the spectral initialization machinery of
[`InternalModesSpectral`](/internal-modes/classes/numerical-solvers/internalmodesspectral)
to build smooth representations of `\rho(z)` and `N^2(z)`, then
applies the WKB asymptotic formulas for the fixed-`\omega`
eigenvalue problem.

This class is most useful for asymptotic comparison, for exploring
turning-point structure, and for the analytical Airy-style
approximations discussed around Sections 2.3, 4.3, and 4.4 of
Early, Lelong, and Smith (2020).

```matlab
im = InternalModesWKB(rho=rho, zIn=zIn, zOut=zOut, latitude=latitude);
[F, G, h, k] = im.ModesAtFrequency(5*im.f0);
```




## Topics
+ Create and initialize modes
  + [`InternalModesWKB`](/internal-modes/classes/analytical-models/internalmodeswkb/internalmodeswkb.html) Initialize the WKB approximation solver.
+ Inspect analytical solutions
  + [`ModesAtFrequencyAiry`](/internal-modes/classes/analytical-models/internalmodeswkb/modesatfrequencyairy.html) Return the turning-point-aware Airy approximation for fixed `\omega`.
  + [`N2_zLobatto`](/internal-modes/classes/analytical-models/internalmodeswkb/n2_zlobatto.html) Return the spectral Lobatto samples of `N2(z)`.
  + [`zLobatto`](/internal-modes/classes/analytical-models/internalmodeswkb/zlobatto.html) Return the spectral Lobatto grid in physical depth coordinates.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`Diff1_zCheb`](/internal-modes/classes/analytical-models/internalmodeswkb/diff1_zcheb.html) Return the spectral first-derivative operator in Chebyshev space.
  + [`ModesAtFrequencyApproximatedAiry`](/internal-modes/classes/analytical-models/internalmodeswkb/modesatfrequencyapproximatedairy.html) Return the simplified WKB Airy approximation for fixed `\omega`.
  + [`SurfaceModesAtWavenumberAlt`](/internal-modes/classes/analytical-models/internalmodeswkb/surfacemodesatwavenumberalt.html) Return the alternate WKB approximation for the surface SQG mode.


---