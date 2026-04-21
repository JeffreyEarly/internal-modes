---
layout: default
title: Class documentation
nav_order: 3
has_children: true
permalink: /classes
mathjax: true
---

# Class Documentation

This section is the API reference for the `internal-modes` class hierarchy. Use it when you want the exact constructor, property, method, or governing equation for one of the solver classes.

The class pages are organized so shared behavior lives on the highest class that defines it. In particular, [`InternalModesBase`](./core-classes/internalmodesbase) documents the common mode-normalization and boundary-condition contract, while solver subclasses focus on what changes in the numerical or analytical formulation.

## Choose a starting point

| Class group | Use it when |
| --- | --- |
| [`Core classes`](./core-classes) | you want the wrapper API or the shared base-class contract |
| [`Numerical solvers`](./numerical-solvers) | you want direct control over the spectral, adaptive-grid, density-coordinate, WKB-coordinate, or finite-difference solvers |
| [`Analytical and asymptotic models`](./analytical-models) | you want closed-form constant/exponential solutions or WKB approximations |
| [`Supporting types`](./supporting-types) | you want the normalization and boundary-condition enumerations |

## Shared notation

The reference pages follow the notation in Early, Lelong, and Smith (2020):

| Symbol | Meaning | API name |
| --- | --- | --- |
| $$F_j(z)$$ | horizontal-velocity vertical structure function | `F` |
| $$G_j(z)$$ | vertical-velocity or density vertical structure function | `G` |
| $$h_j$$ | equivalent depth or eigendepth | `h` |
| $$N^2(z)$$ | buoyancy frequency squared | `N2` |
| $$f_0$$ | Coriolis parameter at the chosen latitude | `f0` |
| $$K$$ | horizontal wavenumber magnitude | `k` |
| $$\omega$$ | wave frequency | `omega` |
| $$z \in [-D, 0]$$ | physical depth coordinate | `z`, `zIn`, `zOut` |
| $$s(z)$$ | stretched coordinate used by some solvers | `x_function`, `xLobatto` |

For the main vertical eigenvalue problems, the package follows the manuscript's two practical forms for the `G` modes:

$$
\partial_{zz} G_j = -\frac{N^2 - \omega^2}{g h_j} G_j
$$

for fixed $$\omega$$, and

$$
\partial_{zz} G_j - K^2 G_j = -\frac{N^2 - f_0^2}{g h_j} G_j
$$

for fixed $$K$$.

## Reading the reference

- [`InternalModes`](./core-classes/internalmodes) is the usual entry point when you want one API that chooses among the concrete implementations.
- [`InternalModesSpectral`](./numerical-solvers/internalmodesspectral) is the reference page for the shared spectral machinery used by the stretched-coordinate subclasses.
- [`InternalModesDensitySpectral`](./numerical-solvers/internalmodesdensityspectral) and [`InternalModesWKBSpectral`](./numerical-solvers/internalmodeswkbspectral) document the two main stretched coordinates from Sections 4.2 and 4.3 of the manuscript.
- [`InternalModesAdaptiveSpectral`](./numerical-solvers/internalmodesadaptivespectral) documents the adaptive multi-region frequency-grid strategy from Section 4.4.
- [`InternalModesConstantStratification`](./analytical-models/internalmodesconstantstratification) and [`InternalModesExponentialStratification`](./analytical-models/internalmodesexponentialstratification) are the main analytical benchmarks used throughout the package and its smoke tests.
