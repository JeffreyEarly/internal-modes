---
layout: default
title: InternalModesWKBHydrostatic
has_children: false
has_toc: false
mathjax: true
parent: Analytical and asymptotic models
grand_parent: Class documentation
nav_order: 4
---

#  InternalModesWKBHydrostatic

Compute hydrostatic WKB mode approximations from a spectrally resolved stratification.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef InternalModesWKBHydrostatic < InternalModesSpectral</code></pre></div></div>

## Overview

`InternalModesWKBHydrostatic` uses the spectral initialization path of
[`InternalModesSpectral`](/internal-modes/classes/numerical-solvers/internalmodesspectral)
to resolve $$N^2(z)$$, then applies a hydrostatic WKB approximation for
the fixed-$$\omega$$ problem. This class is useful as an asymptotic
comparison tool rather than as the primary production solver.

In this approximation, the vertical phase coordinate is built from
the positive part of $$N(z) - \omega$$, and the modal depth is
estimated from

$$
h_j = \frac{1}{g}\left(\frac{d}{j\pi}\right)^2,
$$

where $$d$$ is the accumulated hydrostatic WKB phase over the
oscillatory region.

```matlab
im = InternalModesWKBHydrostatic(rho=rho, zIn=zIn, zOut=zOut, latitude=latitude);
[F, G, h, k] = im.ModesAtFrequency(5*im.f0);
```




## Topics
+ Create and initialize modes
  + [`InternalModesWKBHydrostatic`](/internal-modes/classes/analytical-models/internalmodeswkbhydrostatic/internalmodeswkbhydrostatic.html) Initialize the hydrostatic WKB approximation.


---