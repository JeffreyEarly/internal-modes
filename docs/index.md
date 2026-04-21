---
layout: default
title: Home
nav_order: 1
description: "Compute vertical modes and SQG modes for arbitrary stratification"
permalink: /
---

# Internal Modes

`Internal Modes` provides MATLAB tools for computing baroclinic vertical modes and SQG boundary modes for arbitrary stratification.

The package includes spectral, adaptive-grid, finite-difference, analytical, and WKB-based solvers, all organized around the vertical eigenvalue problems described in Early, Lelong, and Smith (2020).

If you use these classes in published work, please cite:

- J. Early, M. P. Lelong, and K. S. Smith, *Fast and Accurate Computation of Vertical Modes*, *Journal of Advances in Modeling Earth Systems* (2020), [doi:10.1029/2019MS001939](https://doi.org/10.1029/2019MS001939).

## Start here

- [Installation](installation)
- [Class documentation](classes)
- [Version history](version-history)

## Main entry points

| Class | Use it for |
| --- | --- |
| [`InternalModes`](classes/core-classes/internalmodes) | wrapper-based initialization that chooses an implementation and exposes the standard workflow |
| [`InternalModesSpectral`](classes/numerical-solvers/internalmodesspectral) | direct spectral solution of the depth-coordinate EVP |
| [`InternalModesWKBSpectral`](classes/numerical-solvers/internalmodeswkbspectral) | WKB stretched-coordinate spectral solution |
| [`InternalModesAdaptiveSpectral`](classes/numerical-solvers/internalmodesadaptivespectral) | adaptive multi-region solution of the frequency-constant EVP |
| [`InternalModesConstantStratification`](classes/analytical-models/internalmodesconstantstratification) | closed-form constant-stratification solutions |
| [`InternalModesExponentialStratification`](classes/analytical-models/internalmodesexponentialstratification) | closed-form exponential-stratification solutions |
