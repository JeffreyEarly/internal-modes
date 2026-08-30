---
layout: default
title: Analytical bases
parent: Class documentation V2
nav_order: 3
has_children: true
permalink: /classes-v2/analytical-bases
mathjax: true
---

Reference pages for analytical solution families, exact internal-mode bases, and SQG boundary modes under constant and exponential stratification profiles.

The constant and exponential analytical families support the public generalized-energy geostrophic APV descriptor `IMInternalModes.geostrophicAPVModes`. Both apply either the free-surface or rigid-lid convention, handle signed finite, zero (Dirichlet), and positive-infinite endpoint accelerations, and return modes in ascending eigenvalue order $$1/h$$: up to two negative modes, an optional exact zero mode with `h=Inf`, then positive modes.

For $$N^2=N_0^2$$, the exact catalog uses hyperbolic negative branches, an affine zero branch, and trigonometric positive branches. For $$N^2(z)=N_0^2\exp(2z/b)$$, it uses the corresponding exact modified-Bessel, integrated-zero, and ordinary-Bessel branches. APV bases use the volume-only `depth` normalization by default. Roots, normalized determinant residuals, branch classifications, endpoint inertia, and the direct `g0`, `gd`, and `surfaceBoundary` contract are exposed in basis metadata and mode-selection diagnostics.
