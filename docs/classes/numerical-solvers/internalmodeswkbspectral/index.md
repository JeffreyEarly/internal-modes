---
layout: default
title: InternalModesWKBSpectral
has_children: false
has_toc: false
mathjax: true
parent: Numerical solvers
grand_parent: Class documentation
nav_order: 2
---

#  InternalModesWKBSpectral

Solve the vertical EVP on a WKB stretched coordinate with Chebyshev collocation.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef InternalModesWKBSpectral < InternalModesSpectral</code></pre></div></div>

## Overview

`InternalModesWKBSpectral` implements the WKB-coordinate spectral
method described in Section 4.3 of Early, Lelong, and Smith (2020).
It introduces the stretched coordinate

$$
s(z) = \int_z^D N(z') \, dz',
$$

and solves the transformed fixed-$$K$$ and fixed-$$\omega$$
eigenproblems in $$s$$, with

$$
F_j = h_j N \, \partial_s G_j.
$$

Compared with `InternalModesSpectral`, this class concentrates grid
resolution where stratification is strong while preserving the
public constructor contract used by downstream packages.

```matlab
im = InternalModesWKBSpectral(rho=rho, zIn=zIn, zOut=zOut, latitude=latitude, nEVP=257);
[F, G, h, omega] = im.ModesAtWavenumber(2*pi/1000);
```




## Topics
+ Create and initialize modes
  + [`InternalModesWKBSpectral`](/internal-modes/classes/numerical-solvers/internalmodeswkbspectral/internalmodeswkbspectral.html) Initialize the WKB-coordinate spectral solver.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`Nz_function`](/internal-modes/classes/numerical-solvers/internalmodeswkbspectral/nz_function.html) Derivative of $$N(z)$$ used when assembling the stretched-coordinate EVP.
  + [`Nz_xLobatto`](/internal-modes/classes/numerical-solvers/internalmodeswkbspectral/nz_xlobatto.html) $$\partial_z N$$ sampled on the Lobatto grid in the WKB coordinate.


---