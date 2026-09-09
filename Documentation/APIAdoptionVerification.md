# Simplified API integration verification

The provider integration retains three new public classes: `IMBasisCollection`, `IMProjection`, and `IMProductAssessmentPlan`. Projection descriptors are hidden implementation data; assessments are ordinary structures with separate checking functions. Released transform APIs remain available. The final runtime source was identical before and after reconstructing the three review branches.

## Numerical checks

The complete `UnitTestsV2` suite passed **348 of 348 tests**, with no failures or incomplete tests. This includes the scientific corrections that were absent from the initially submitted API branches: rotating-wave normalization, consistent stretched-coordinate calculus, row equilibration and finite-mode filtering, and localized analytical APV roots.

The foundation passed 32 focused tests and its executable example. The independently restacked foundation branch also passed all 32 tests. The restacked bulk branch passed its 13 construction/evaluation tests. Product assessment passed 20 focused tests and its cosine-product example, including the frozen mixed-product numerical controls. Coverage includes unequal collection-page counts, normalization snapshots, provenance, endpoint requirements, signed physical duals, positive norms, reference qualification, structural zeros, and rejection before expensive callbacks when reservations exceed the budget.

Reproduce the complete numerical gate with the provider and its manifest dependencies on a clean MATLAB path:

```matlab
results = runtests(fullfile(providerRoot,"UnitTestsV2"));
assertSuccess(results);
```

Generated class documentation was rebuilt separately at each review branch's API level. InternalModes does not currently define a `docs:check` build task; verification uses its documentation builder, generated API-name inspection, and whitespace/source-scope checks. The standalone checking functions have MATLAB help and are described in [AssessmentAPI.md](AssessmentAPI.md) and [ProductInventoryAPI.md](ProductInventoryAPI.md).

## WVM construction evidence

The initial consumer migration passed 22 focused construction, transform, linear-evolution, and resolution-transfer tests, including physical source projection, independent family counts, energy, output, and restart. It fills the existing WVM arrays from exact bulk solves; no provider object graph is persisted.

A bounded benchmark used 33 already-distinct wavenumbers including zero, `nEVP=64`, four wave modes, three inertial modes, and 65 evaluation depths. After one warmup and three measurements, the median ratio of scalar to bulk solve/evaluation time was 2.3902. Both paths performed 33 exact solves; the bulk path used one shared preparation. This is a case-specific construction measurement, excludes balanced-family construction, and is not a general speedup guarantee.

## Existing physical-residual limitation

The corrected provider baseline already fails one WVM long-external-wave vertical-momentum check: approximately `2.41209e-7` against `2e-7` for exponential stratification, horizontal wavelength 100 km, depth 1000 m, `nEVP=64`, and 129 evaluation points. No tolerance was changed during adoption.

A bounded baseline diagnostic varied `nEVP` from 48 to 128 and evaluation counts from 65 to 257. At 129 evaluation points, the long external mode's residual grew from `9.89e-8` to `3.12e-6`; evaluation refinement from 129 to 257 changed it little. Internal modes in the same case had residuals around `1e-11` or smaller. This is consistent with a native-solve/pressure-differentiation conditioning or roundoff limitation; it does not establish a new physical correction. [InternalModes #10](https://github.com/JeffreyEarly/internal-modes/issues/10) retains ownership of solve-quality diagnostics. Passing projection/product tests does not claim complete physical qualification.

Full per-wavenumber WVM retained counts remain the high-priority follow-up [WVM #414](https://github.com/JeffreyEarly/wave-vortex-model/issues/414). The collection layer preserves unequal counts now. Nearby-wavenumber approximation and thermal-mode research are outside this adoption.
