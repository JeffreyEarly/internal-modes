---
layout: default
title: Version History
nav_order: 100
---

# Version History

## [1.2.0] - 2026-04-20
- modernized the non-spectral constructors to explicit name-value initialization while preserving the `InternalModesSpectral` and `InternalModesWKBSpectral` constructor contracts
- updated the `InternalModes` wrapper to forward constructor options through the new named-argument paths without breaking wrapper-level properties such as diagnostics and boundary-condition settings
- renamed the former script-only `UnitTests` collection to `Examples`, repaired the example scripts and local asset loading, and added a real `matlab.unittest` smoke suite under `UnitTests`
- refreshed direct constructor call sites and documentation snippets to match the new constructor forms and current spline dependency APIs

## [1.1.0] - 2026-04-09
- updated spline-based interpolant construction to the `Distributions` 2.0 named-argument API
- raised the `SplineCore` dependency floor to `^2.0` and added a direct `Distributions ^2.0` dependency for the package's direct `NormalDistribution` usage

## [1.0.1] - 2025-12-09
- Initial CI release
