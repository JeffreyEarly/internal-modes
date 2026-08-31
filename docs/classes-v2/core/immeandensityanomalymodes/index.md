---
layout: default
title: IMMeanDensityAnomalyModes
has_children: false
has_toc: false
mathjax: true
parent: Core
grand_parent: Class documentation V2
nav_order: 7
---

#  IMMeanDensityAnomalyModes

Describe generalized-energy mean-density-anomaly modes.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMMeanDensityAnomalyModes < IMInternalModes</code></pre></div></div>

## Overview

`IMMeanDensityAnomalyModes` describes the `G`-form problem

$$
-G_j''(z)=\frac{N^2(z)}{g h_j}G_j(z).
$$

A finite surface acceleration imposes
$$g h_jG_j'(z_s)=g_0G_j(z_s)$$, and a finite bottom
acceleration imposes
$$g h_jG_j'(z_b)=-g_dG_j(z_b).$$ Positive infinity makes the
corresponding endpoint inactive and imposes Dirichlet data. Finite
zero therefore gives the Neumann limit.

The solved `G` modes use the signed generalized-energy product

$$
\langle G_i,G_j\rangle=
\frac{1}{g}\int_{z_b}^{z_s}N^2G_iG_j\,dz
+\frac{g_0}{g}G_i(z_s)G_j(z_s)
+\frac{g_d}{g}G_i(z_b)G_j(z_b),
$$

with inactive endpoint terms omitted. The companion pressure mode is
the surface-referenced diagnostic

$$
F_j(z)=\frac{1}{g}\int_z^{z_s}N^2(z')G_j(z')\,dz',
\qquad F_j(z_s)=0.
$$

Construct this descriptor through
`IMInternalModes.meanDensityAnomalyModes`.




## Topics
+ Create internal-mode EVPs
  + [`IMMeanDensityAnomalyModes`](/internal-modes/classes-v2/core/immeandensityanomalymodes/immeandensityanomalymodes.html) Create a generalized-energy mean-density-anomaly EVP.
+ Inspect internal-mode configuration
  + [`activeEndpoints`](/internal-modes/classes-v2/core/immeandensityanomalymodes/activeendpoints.html) Active endpoint identities in canonical surface-bottom order.
  + [`g0`](/internal-modes/classes-v2/core/immeandensityanomalymodes/g0.html) Surface generalized-energy acceleration.
  + [`gd`](/internal-modes/classes-v2/core/immeandensityanomalymodes/gd.html) Bottom generalized-energy acceleration.
+ Inspect internal-mode inner products
  + [`innerProduct`](/internal-modes/classes-v2/core/immeandensityanomalymodes/innerproduct.html) Return the mean-density-anomaly `F` or `G` inner product.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`makeBasisSet`](/internal-modes/classes-v2/core/immeandensityanomalymodes/makebasisset.html) Create a mean-density-anomaly basis set.


---