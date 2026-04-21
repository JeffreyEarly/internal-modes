---
layout: default
title: InternalModesDensitySpectral
has_children: false
has_toc: false
mathjax: true
parent: Numerical solvers
grand_parent: Class documentation
nav_order: 4
---

#  InternalModesDensitySpectral

Solve the vertical EVP on a density-stretched coordinate with Chebyshev collocation.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef InternalModesDensitySpectral < InternalModesSpectral</code></pre></div></div>

## Overview

`InternalModesDensitySpectral` implements the density-coordinate
formulation discussed in Section 4.2 of Early, Lelong, and Smith
(2020). It uses the stretched coordinate

$$
s(z) = -\frac{g}{\rho_0}\rho(z) + g,
$$

so that

$$
F_j = h_j N^2 \partial_s G_j.
$$

This coordinate concentrates points where the density varies most
rapidly and therefore requires a monotonic background density
profile.

```matlab
im = InternalModesDensitySpectral(rho=rho, zIn=zIn, zOut=zOut, latitude=latitude, nEVP=257);
[F, G, h, omega] = im.modesAtWavenumber(2*pi/1000);
```




## Topics
+ Create and initialize modes
  + [`InternalModesDensitySpectral`](/internal-modes/classes/numerical-solvers/internalmodesdensityspectral/internalmodesdensityspectral.html) Initialize the density-coordinate spectral solver.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`N2z_xLobatto`](/internal-modes/classes/numerical-solvers/internalmodesdensityspectral/n2z_xlobatto.html) $$\partial_z N^2$$ sampled on the Lobatto grid in the density coordinate.


---