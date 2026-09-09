# Resolved-mode convergence and counts by wavenumber

`assessModeConvergence(candidate,reference,z,weights)` compares the same physical modes computed at two explicitly chosen eigenproblem resolutions. It returns ordinary structs and measurement tables. It performs no eigensolve, chooses no tolerance or retained count, and changes no scientific modes. A caller can reuse an already prepared reference across assessments. This is the shared comparison used by the WVM advisory workflow; sampled-transform and quadratic-product measurements remain separate operations.

## Prepared inputs and identities

Each input contains `identity`, `values` and `provenance`. `identity` identifies the family, exact scientific `columnLabels`, normalization and `zDomain`; include `kappa` for wave families. All identity fields other than the column labels must agree exactly. `values` contains sampled physical fields, such as `F` and `G`, with one column per label. Optional `derivatives` contains their physical vertical derivatives, and optional `equivalentDepths` and `eigenvalues` contain scalar mode properties. `provenance` records the actual solver, coordinate kind, eigenproblem resolution and source revision, with any additional reference information the caller needs.

The comparison matches exact labels in candidate order. It reports missing or duplicated labels as inconclusive; it does not guess a match from nearby eigenvalues. A common sign aligns each paired mode across all its fields and derivatives, preserving their physical polarization. The inputs must describe the same physical problem on the same increasing physical sample points with positive integration weights.

```matlab
identity = struct(family="waves",columnLabels=string(basis.modeNumber),normalization=string(basis.normalization),zDomain=basis.zDomain,kappa=kappa);
candidate = struct(identity=identity,values=struct(G=basis.G(z)),equivalentDepths=basis.h,provenance=struct(nEVP=64));
reference = struct(identity=identity,values=struct(G=refined.G(z)),equivalentDepths=refined.h,provenance=struct(nEVP=128));
assessment = assessModeConvergence(candidate,reference,z,weights);
disp(assessment.matches)
disp(assessment.measurements)
```

This short illustration assumes that the independently solved `basis` and `refined` have identical labels. When they differ, each input must carry its own actual labels. The runnable example below includes the complete construction, derivative evaluation, provenance and interpretation.

## Measurements and their limits

For a sampled physical field, define the positive discrete norm by $$\|v\|_w^2 = \sum_i w_i |v(z_i)|^2.$$ After common sign alignment, shape and derivative measurements are weighted relative L² discrepancies against the reference. Joint H¹ discrepancies use domain depth $$D = z_s - z_b$$:

$$e_{H^1}(v) = \frac{\sqrt{\|v - v_r\|_w^2 + D^2 \|v_z - v_{r,z}\|_w^2}}{\sqrt{\|v_r\|_w^2 + D^2 \|v_{r,z}\|_w^2}}.$$

This joint norm remains useful when a nearly constant field has a roundoff-sized derivative. The derivative-only relative discrepancy is still reported; it can be large in that situation without implying a comparably large field error. No smallness threshold based on unrelated modes hides an individual mode's error. A zero reference norm gives zero discrepancy only when the difference is also zero; otherwise the measurement is infinite.

Finite equivalent depths are compared using $$e_h = |h - h_r|/|h_r|,$$ with the analogous expression for eigenvalues. Equal zeros and equal signed infinities agree. Finite-versus-infinite values, opposite infinities and NaNs are inconclusive. Optional quantities not supplied are marked `notRequested`.

The result includes per-mode matching/orientation, long-form measurements, input identities/provenance, coverage and comparison cost. A `measured` row means that a comparison was possible. Agreement at two resolutions supplies convergence evidence, not an absolute error bound. The function does not independently qualify the reference solve or supplied quadrature; its coverage states those limits explicitly. Accurate solver regressions, convergence under further refinement and application-specific checks continue to provide complementary evidence.

## Different wave counts at each wavenumber

`IMSolver.solveWaveModesAtWavenumbers` accepts either scalar `nModes` or a row aligned with the positive entries of the requested `kappa` row, including repeats. Repeated identical wavenumbers must request identical counts. Zero-wavenumber pages instead use the explicit independent scalar `nInertialModes`.

```matlab
kappa = [2e-4 0 1e-4 2e-4];
[collection,costs] = solver.solveWaveModesAtWavenumbers(kappa,N2=N2,zDomain=[-1000 0],nModes=[4 7 4],nInertialModes=5);
```

Each distinct wavenumber is solved once, using shared preparation and the exact requested count. The collection preserves the original basis, labels and page mapping. Evaluating pages with different column counts together requires an explicit common column selection, or separate calls; rectangular outputs do not silently pad or truncate the scientific bases. Wave frequency signs and model storage conventions remain the consumer's responsibility.

## Reproduce the figure

With this development version of InternalModes and its declared dependencies on the MATLAB path, run:

```matlab
addpath ExamplesV2
result = waveModeCountsByWavenumber();
% A named directory must not already exist:
result = waveModeCountsByWavenumber(Nz=33,outputDirectory="wave-counts-33");
```

The [example source](../ExamplesV2/waveModeCountsByWavenumber.m) uses a 1 km square periodic domain, 16 × 16 horizontal points, depth 1000 m, and a fixed 25-point WKB-stretched Chebyshev–Lobatto physical quadrature. The profile is $$N^2(z) = N_0^2 \exp(2z/b)$$ with $$N_0 = 0.01\,\mathrm{s}^{-1}$$ and $$b = 400\,\mathrm{m}$$; $$f_0 = 10^{-4}\,\mathrm{s}^{-1}$$ and $$g = 9.81\,\mathrm{m}\,\mathrm{s}^{-2}$$. It solves every distinct positive Fourier radius of that horizontal grid, requesting 24 modes including the external wave, at eigenproblem resolutions 96 and 144. The independently requested inertial inventory contains 12 modes.

The plot separates three contiguous prefixes: modes passing the relative equivalent-depth and both F/G H¹ comparisons at tolerance $$10^{-5}$$; modes passing the sampled G Gram test at tolerance 1%; and modes passing both. The physical G metric retains its $$(N^2-f_0^2)/g$$ interior weight and free-surface endpoint contribution. The zero-wavenumber inertial comparison uses its separate physical F identity $$\int F_i F_j\,dz = h_i\delta_{ij}$$. No density displacement, mode normalization or wave polarization is redefined to improve the plot.

Candidate-ceiling markers distinguish a lower bound from a measured limit. Inconclusive comparisons are marked separately. The default integration run found all 24 candidate waves converged, with the sampled grid supporting between 13 and 8 leading modes over the positive wavenumbers; all 12 requested inertial modes passed. This is an example result with its explicit tolerances, not a universal retention rule. Counts need not vary monotonically because different modal shapes can interact differently with a fixed quadrature.

Each run writes PNG/PDF figures, wave/inertial count CSVs, per-mode convergence and per-prefix Gram measurements, the physical grid and weights, and JSON configuration/source provenance. Without an output-directory argument it chooses a unique temporary directory and prints the path; existing output directories are never overwritten. Increasing candidate count or eigenproblem resolution can test the plotted ceiling, and changing vertical point count tests the fixed-grid limitation. No part of this figure certifies quadratic aliasing or leakage from arbitrary excluded modes.

The [recorded default figure](ModeConvergenceExample/wave-counts.png), [vector PDF](ModeConvergenceExample/wave-counts.pdf), [counts](ModeConvergenceExample/wave-counts.csv) and [provenance](ModeConvergenceExample/provenance.json) were generated from the committed implementation recorded in that provenance. Full per-mode measurements are reproducible by running the example.
