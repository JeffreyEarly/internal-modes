---
layout: default
title: Normalization
has_children: false
has_toc: false
mathjax: true
parent: Supporting types
grand_parent: Class documentation V2
nav_order: 1
---

#  Normalization

Enumerate the supported modal normalization conventions.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef Normalization</code></pre></div></div>

## Overview

`Normalization` collects the normalization choices used throughout
the `internal-modes` hierarchy. Following Section 2.4 of Early,
Lelong, and Smith (2020), the main wave-mode norms are the
`unity` norm based on the active EVP inner product, the `kConstant`
norm based on the fixed-$$K$$ orthogonality relation, and the
`omegaConstant` norm based on the fixed-$$\omega$$ relation.

The valid values are:

- `Normalization.unity` for unit norm under the active EVP inner product
- `Normalization.kConstant` for the manuscript's $$K$$-constant norm
- `Normalization.omegaConstant` for the manuscript's
  $$\omega$$-constant norm
- `Normalization.uMax` to scale each mode so $$\max F_j = 1$$
- `Normalization.wMax` to scale each mode so $$\max G_j = 1$$
- `Normalization.surfacePressure` to scale by the raw surface value
  of $$F$$ so $$F_j^\mathrm{surfacePressure}(z_\mathrm{surface})=1$$
- `Normalization.geostrophic` for the near-geostrophic interior mode
  normalization used by some helper workflows

```matlab
im = InternalModes(rho, zIn, zOut, latitude);
im.normalization = Normalization.kConstant;
```




## Topics


---