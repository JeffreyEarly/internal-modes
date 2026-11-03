# Prescribed vertical sampling verification

This increment follows the PR #66 integration baseline `43d1e4534ab5ffd3c68f6e91da9d5fba00115954`. It adds vertical point/cell sampling and a portable diagnostic-layout control for [#60](https://github.com/JeffreyEarly/internal-modes/issues/60). It does not add WVM physical projectors or change the shared recovery engine.

## Portable controls

`UnitTestsV2/IMVerticalSamplingTests.m` contains nine methods. Its committed CSV supplies six synthetic, staggered F-point/G-cell-mean observations on [-2, 0] m, in an arbitrary fixed order, with no endpoint samples or WVM dependency. Independent Fourier primitives supply expected point, mean and integral values. A selected harmonic 137 tests high-frequency integration and explicit budget failure.

Other controls compare neutral-label cells with clipped sine integrals and classical piecewise derivatives, and finite-difference F/G cells with exact integration of the independently reconstructed PCHIP coefficients on their distinct knot sets. Selection, normalization, dense zero-APV rotations and repeated wavenumber pages preserve the evaluated matrix and scientific identities within numerical tolerance.

The shared recovery fixture keeps every scientific coefficient active. Its constant F/zero G sector is recoverable from the staggered F/G layout; G-only observations fail explicitly. Complex coefficients recover in the declared finite span. A separate nullspace construction produces an off-span field with identical observations but a different value at a hidden point, demonstrating the limit of observation-rank evidence.

Snapshot controls preserve row/column identities and prepared matrices after subsequent basis changes, and count no additional EVP solves during sampling or repeated transform application. The executable `ExamplesV2/PrescribedVerticalSampling.m` additionally exercises a 5 × 2 × 3 complex coefficient batch with ten irregular/staggered observations; its maximum coefficient round-trip error was 4.45e-16. This measures prepared finite-span recovery, not independent physical projection.

## Independent mathematical and numerical audit

A separate reviewer used references beyond the committed regression fixtures. Errors below are absolute control-specific matrix infinity norms; their units and normalization differ across rows, so they are not one common physical error measure.

| Control | Measured error or result |
| --- | --- |
| Selected DD harmonics 137/160, partial/support-crossing G integrals | 1.76e-15 |
| Same DD columns, APV cell means | 1.84e-12 |
| Same DD columns, classical Gzz integral | 6.94e-10 |
| Fixed 32-point Gauss comparison on a high-harmonic cell | 0.313 error, illustrating under-resolution |
| Actual finite-difference companion F/G knots, exact piecewise-cubic reference | 2.78e-15 |
| Stored FD Gz integral versus endpoint difference | 0.361 difference, confirming that the literal stored derivative must be integrated |
| Spectral physical-z columns 137/150, polynomial-exact reference | 3.18e-12 |
| Mapped spectral source, independent refined native-coordinate reference | 1.13e-15 |
| Mixed point/cell zero-APV coordinate rotation | 3.33e-16; a 1e6 coordinate-scale ratio also passed |

The audit found a floating-point range defect: on depth [0, 1e-306], a constant F response of 1e-16 had a representable integral near 1e-322, but individual weighted contributions rounded to zero and falsely reported convergence. Preparation now rejects complete positive-mass loss with `IMVerticalSampling:UnrepresentableIntegral`. It also rejects zero/nonfinite positive physical quadrature weights. An independent re-test confirmed the rejection and preserved the corresponding representable mean, 9.999999999999997e-17. A dedicated regression covers this case. The guard detects complete mass loss; it does not certify every subnormal computation.

Mathematics, API ownership, example and test reviews found no remaining required finding after this correction. The sampler owns vertical functionals and their evidence, the continuous basis owns profiles and coordinates, the shared engine owns numerical recovery, and WVM owns its complete physical state/source combinations.

## Verification and limits

The focused run passed all 35 selected methods and the executable example. After the bounded range correction, its new regression passed independently. The frozen-source full V2 run then passed all 709 methods, with zero failures or incomplete tests. Static analysis identified implicit sharing of two loop variables between the point evaluator and a nested integration helper; the helper's local names were isolated afterward. The exact four-line binding-only delta was checked against the frozen source; all other source/fixture hashes remained unchanged. All nine sampling methods and the executable example then passed again. Code Analyzer has zero remaining findings across the ten changed MATLAB files.

Temporary `ClassDocumentation` extraction for the new sampling class and affected common basis succeeded. Markdown local links and fences, whitespace, unchanged package metadata, repository scope and generated-artifact checks passed. This repository has no configured `docs:check` task; no website output or released package snapshot was edited. All MATLAB gates ran in fresh local processes with one computational thread.

To replay the portable suite with InternalModes and its declared dependencies on the path:

```matlab
results = runtests("UnitTestsV2/IMVerticalSamplingTests.m");
assert(all([results.Passed]));
addpath ExamplesV2
result = PrescribedVerticalSampling();
```

Cell integration records successive-rule estimates and the finer rule's absolute-profile mass estimate. Arbitrary unresolved intrinsic features can evade those estimates. Mode-solve accuracy, continuous normalization references, physical quadrature, full WVM projectors, nonlinear products and arbitrary excluded modal content remain separate qualification tasks. In particular, this increment does not replace or qualify the existing analytical-source trapezoidal physical-reference construction described in [common scalar projection](CommonScalarProjection.md).
