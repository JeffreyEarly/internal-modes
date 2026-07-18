---
layout: default
title: UpperBoundary
has_children: false
has_toc: false
mathjax: true
parent: Supporting types
grand_parent: Class documentation
nav_order: 2
---

#  UpperBoundary

Enumerate the supported upper boundary conditions.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef UpperBoundary</code></pre></div></div>

## Overview

pperBoundary` collects the upper-surface boundary-condition labels
ed by the `internal-modes` solvers. The main user-facing values are
igidLid` and `freeSurface`, matching the boundary conditions
scribed in Section 2.3 of Early, Lelong, and Smith (2020).

e valid values are:

`UpperBoundary.rigidLid` for $$G(0) = 0$$
`UpperBoundary.freeSurface` for $$h \partial_z G(0) = G(0)$$
`UpperBoundary.mda`, `UpperBoundary.buoyancyAnomaly`,
`UpperBoundary.geostrophicFreeSurface`, and `UpperBoundary.custom`
for specialized internal workflows
`UpperBoundary.none` when no explicit upper condition should be
imposed

`matlab
 = InternalModes(rho, zIn, zOut, latitude);
.upperBoundary = UpperBoundary.freeSurface;
`




## Topics


---