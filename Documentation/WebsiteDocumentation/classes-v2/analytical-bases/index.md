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

The exponential analytical family also supports the generalized-energy geostrophic APV descriptor named `geostrophicAPVModes`. For $$N^2(z)=N_0^2\exp(2z/b)$$ it solves the positive, zero, and negative equivalent-depth branches exactly, applies either the free-surface or rigid-lid endpoint convention, and handles finite, zero (Dirichlet), and positive-infinite endpoint accelerations. Modes are returned in ascending eigenvalue order $$1/h$$: up to two negative modes, an optional exact zero mode with `h=Inf`, and then the positive modes. The endpoint-inertia calculation and branch classifications are exposed in the basis metadata and mode-selection diagnostics.
