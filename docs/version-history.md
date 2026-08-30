---
layout: default
title: Version History
nav_order: 100
---

# Version History

## Unreleased
- replaced the projected surface-geostrophic API with coefficient-independent `IMGeostrophicZeroAPVModes` boundary-response problems, page-shaped canonical `F/G` bases, explicit physical-energy and endpoint-buoyancy forms, and boundary-depth, surface-buoyancy, and custom symmetric-pencil rotations
- added the public v2 `IMInternalModes.geostrophicAPVModes` generalized-energy APV factory, signed endpoint-limit mappings and metadata, exact numerical zero retention, volume-only `depth` normalization, and exact constant-stratification trigonometric/affine/hyperbolic branches with independent factory, endpoint, normalization, analytical-root, and refined-spectral regression coverage
- added exact exponential-stratification generalized-energy geostrophic APV modes, including positive, zero, and negative equivalent-depth branches, free-surface and rigid-lid conventions, all finite/zero/infinite endpoint limits, signed-spectrum diagnostics, and analytical-versus-numerical regression coverage
- added v2 scalar discrete Galerkin transforms on caller-supplied points and quadrature weights, including forward and inverse matrices, forward and back transform methods, sampled metrics, continuous Gram targets, and transform-quality diagnostics
- added fixed-point v2 quadrature fitting with normalized Gram objectives, nonnegative full-depth constraints, custom least-squares systems, and fitted-versus-geometric diagnostics
- made fixed-point fitting depth-scale robust and added v2 walkthroughs for default and custom scalar discrete-transform workflows
- added v2 `pointsFromModeRoots` grids with automatic auxiliary-mode solves and native Chebyshev root finding for all spectral coordinate choices
- compressed the normalized Gram Frobenius fit to independent symmetric mode pairs and made `discreteTransform` optionally return the associated quadrature-fit provenance when it chooses weights

## [2.0.0] - 2026-05-19
- added the parallel v2 `IMEigenvalueProblem` architecture with physical-coordinate EVP descriptors, structured operators, first-class boundary conditions, coordinate-aware spectral/WKB/density/finite-difference solvers, and `IMBasisSet`
- moved v2 modal normalization ownership to `IMBasisSet` so solved native modes can be reinterpreted with different normalizations without re-solving the EVP
- added v2 analytical solution families for constant and exponential stratification
- renamed v2 wave EVPs to `waveModesAtWavenumber` and `waveModesAtFrequency`, with `hydrostaticGModes` now using the fixed-frequency wave EVP at `omega=0`
- added native-coordinate v2 inner products for spectral, WKB spectral, and density spectral solvers, including endpoint surface and bottom weights outside the interior integral
- added separate `UnitTestsV2` regression coverage and `ExamplesV2` scripts for the new IM stack
- added index policies for active-boundary mode counts, including the manuscript positive-boundary and negative-boundary PE conventions
- added v2 regression tests for coordinate pullbacks, index policies, finite-difference grid ownership, basis-owned normalization, partial-depth Gram matrices, and unsupported analytical-basis operations
- expanded the partial-depth observational projection example with configurable sampling scenarios, Gaussian ensemble coefficient realizations, physical G-mode potential-energy spectra, observed-space completeness diagnostics, and spectral-window plots
- added a hydrostatic quadrature accuracy comparison example that separates mode-shape, root, solver, retained-band, and positive-weight effects in Parseval diagnostics
- added geostrophic normalization support and diagnostics for the analytical constant and exponential stratification solvers, with smoke-test coverage
- updated the physical sampling quadrature optimization example to use geostrophic normalization

## [1.4.0] - 2026-05-08
- added component-role-aware `InternalModesBasis`, `InternalModesTransform`, and `InternalModesProjection` classes for vertical transforms, fixed-grid projections, observational projections, vertical spectra, and annotated persistence
- added vertical-transform examples for wave-vortex-model integration sketches and component-role-aware modal spectra
- added a direct `ClassAnnotations ^1.2.1` dependency for NetCDF-backed annotated persistence through `CAAnnotatedClass`

## [1.3.0] - 2026-05-06
- ported the spectral speedup refactor into the `@InternalModesSpectral` class-folder layout
- added `modesAtQuadraturePoints` for computing resolved modes on mode-adapted quadrature points
- refreshed spectral internal state when `nEVP` changes, including WKB stretched-grid derivative state
- fixed `BSpline` initialization in `InternalModesBase`
- added speedup exploration helpers under `Extras/SpeedupExploration`
- expanded automated smoke coverage for quadrature points, `nEVP` refresh behavior, legacy quadrature wrappers, and spline initialization

## [1.2.0] - 2026-04-20
- modernized the non-spectral constructors to explicit name-value initialization while preserving the `InternalModesSpectral` and `InternalModesWKBSpectral` constructor contracts
- updated the `InternalModes` wrapper to forward constructor options through the new named-argument paths without breaking wrapper-level properties such as diagnostics and boundary-condition settings
- renamed the former script-only `UnitTests` collection to `Examples`, repaired the example scripts and local asset loading, and added a real `matlab.unittest` smoke suite under `UnitTests`
- refreshed direct constructor call sites and documentation snippets to match the new constructor forms and current spline dependency APIs
- modernized the core user-facing mode API to lowerCamel method names such as `modesAtFrequency`, `modesAtWavenumber`, `surfaceModesAtWavenumber`, and `showLowestModesAtFrequency`, while preserving the legacy UpperCamel names as hidden compatibility aliases and updating the docs, examples, and smoke tests accordingly
- modernized the direct eigenmatrix helper API to lowerCamel `eigenmatricesForFrequency` and `eigenmatricesForWavenumber`, while preserving the legacy UpperCamel names as hidden compatibility aliases and updating the developer-facing docs and smoke tests accordingly

## [1.1.0] - 2026-04-09
- updated spline-based interpolant construction to the `Distributions` 2.0 named-argument API
- raised the `SplineCore` dependency floor to `^2.0` and added a direct `Distributions ^2.0` dependency for the package's direct `NormalDistribution` usage

## [1.0.1] - 2025-12-09
- Initial CI release
