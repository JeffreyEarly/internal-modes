# Shared recovery verification

Historical verification of the incoming implementation before the capable-basis integration. Code examples and class names below describe that historical API; current usage is documented in [the integration plan](APISimplificationPlan.md). These results do not replace verification of the integrated implementation.

Evidence for [#57](https://github.com/JeffreyEarly/internal-modes/issues/57), measured on 2026-11-02 against base revision `8717aca878fff72b5e2c19f962e7dc6c8c435888`. The implementation and portable controls accompany this note on the same branch. This is provider recovery qualification, not a new physical-grid or maximum-band claim.

## Environment and reproduction

MATLAB `26.2.0.3433840 (R2026b) Update 1`, `MACA64`, one numerical thread. The authoring repositories for Chebfun, SplineCore, Distributions and ClassAnnotations were on the test path. Local MATLAB commands ran outside the macOS sandbox because of the workspace's Qt CPU-probe restriction.

With those dependencies available and this repository on the path:

```matlab
maxNumCompThreads(1);
addpath('UnitTestsV2','UnitTestsV2/TestSupport','tools');
results = runtests('UnitTestsV2');
summary = validateSharedRecovery();
[measurements,environment] = measureSharedRecovery();
```

The new unit-test classes are `IMRecoveryKernelTests`, `IMObservationTransformTests` and `IMInternalModesRecoveryIntegrationTests`. Existing projection, scalar-policy, quadrature, assessment, source and root tests exercise the migrated consumers. The controls construct their matrices in source and need no local manuscript/archive fixture.

## Independent audits

Separate reviewers audited production code they did not author:

- Numerical/mathematical review checked scaled QR least squares, dense complex C/D/F covariance, row-space preparation of weak equations, principal-angle diagnostics, activity and extended-span selection. All 21 kernel methods passed independently. Forty deterministic complex reference controls had maximum relative error `2.11e-13`; thirty dense changes of physical test basis had maximum map change `2.51e-15`. The portable `validateSharedRecovery` tool preserves these reference calculations.
- Framework review checked generic/scalar/aligned class ownership, physical source semantics, cached preparation/application, prefix/reweight routing, caller metadata and composed output views. Independent adapter tests and nested/empty/failed-view controls passed. No material framework finding remained.
- Physical-projection review independently compared dense signed normalization, unit changes, complex prescribed systems, prefix reference forms, physical-only operation, semidefinite fits, inactive/empty sectors and physical/recovery separation. Overflow and complete-column underflow of a normalized physical dual now report unavailable matrices, while independent positive-fit recovery remains usable. The added regression and independent ordinary controls passed.

Auditors exposed and the implementation corrected finite-range normalization defects, an unnecessary observation-square residual allocation and physical-dual availability errors before handoff. The final core audit covered SHA-256 `0b9187ba022866139894157d8debcc71789e9b8facbc50cfa8ee76ce49d7a52a`.

## Verification ledger

- Clean-base focused baseline: 80 of 81 methods passed. `IMDiscreteTransformAssessmentTests/assessmentReproducesFixedRuleAndNestedFit` failed before these changes and passed in the affected integration run.
- Affected integration run: 181 of 183 methods passed initially. The two failures were strict-policy message assertions; restoring accepted count, first failing value and tolerance detail made both targeted reruns pass.
- Later numerical changes: all 12 affected weak/algebraic methods passed, including a 20,000-observation/two-column control and explicit failure outside representable whitening range. The independent final kernel rerun passed all 21 methods.
- Later physical availability change: the expanded overflow/underflow regression, ordinary independent physical controls and both provenance/view tests passed.
- Final frozen-source affected integration: **186 of 186 methods passed**, including all later numerical/physical fixes, scalar/F/G/root transforms, quadrature, source, product and assessment consumers. The portable independent audit also passed, and its final Code Analyzer check reported zero messages.
- Inline API documentation generated 168 Markdown pages for six affected classes into temporary output. No generated site files were edited. The class registry includes the new scalar adapter.
- Code Analyzer covered all changed MATLAB files. The only existing notices were four unchanged-code notices in `IMBasisSet` (class comparison, array growth and two property-name notices); new/changed recovery code was clean. The benchmark's `whos` payload has an explicit unused-variable suppression because its use is by variable name. Files changed after that pass were checked separately.

- Full V2 run: **612 of 614 methods passed**. The two failures were reproduced individually on an export of the exact clean base revision: `IMEigenvalueProblemRefactorTests/developerBasisMethodsArePublicAndCallable` errors because its analytical solution and EVP stratification differ, and `IMSmoothedBandwidthCoordinateTests/constantControlKeepsTheGridAndPhysicalSpectrum` has a `qReference` mismatch of about `1.81e-7` against a `2e-11` relative tolerance. Both source areas and tests are unchanged here. The later weak-memory/range and physical-dual availability changes received the focused and independent reruns recorded above.
- Final whitespace, affected Markdown local links, repository scope, unchanged package/fixture manifests and absence of generated-site/package-snapshot changes were checked. This repository has no configured `docs:check` command; temporary class generation and focused local-link validation were used instead.

## Setup, application and payload costs

The [raw measurements](SharedRecoveryCosts.csv) come from `measureSharedRecovery`, using deterministic complex orthonormal sample frames, full extended-fixture amplitudes spanning `1e-2` to `1e2`, explicitly declared D=I and 256 independent fields per application. Retained matrices use the leading amplitude subset. Each policy is prepared once per application measurement; `timeit` measures setup separately. Extended recovery uses 36 or 144 columns before retaining 24 or 96. No EVP solve, grid design, quadrature fit or fixture generation is included.

| Observations / retained | Policy | Setup (ms) | Forward 256 fields (ms) | Back 256 fields (ms) | Value payload (MiB) |
| --- | --- | ---: | ---: | ---: | ---: |
| 65 / 24 | sampledLeastSquares | 3.28 | 0.064 | 0.059 | 0.132 |
| 65 / 24 | correctedPhysicalDual | 4.70 | 0.031 | 0.032 | 0.160 |
| 65 / 24 | extendedSpan | 6.60 | 0.014 | 0.022 | 0.269 |
| 195 / 96 | sampledLeastSquares | 13.01 | 0.090 | 0.099 | 1.243 |
| 195 / 96 | correctedPhysicalDual | 22.15 | 0.091 | 0.099 | 1.537 |
| 195 / 96 | extendedSpan | 24.10 | 0.091 | 0.100 | 2.599 |

These are illustrative local timings, not a speedup claim or policy ranking. Application performs cached matrix products, and timing variation at these small sizes is substantial. `valuePayloadBytes` is MATLAB `whos` accounting for an explicit struct containing A, R, B, the resolved specification and diagnostics. It excludes object/runtime overhead and is neither unique allocated heap memory nor process RSS. Full matrix factors and extended synthesis provenance contribute to that payload.

## Capability and range limits

Correct sample recovery does not qualify physical quadrature or arbitrary omitted-input suppression. The rounded-Gram control demonstrates observable positive recovery while retaining failed physical Gram evidence. Noise optimality is asserted only for the declared F=D objective. No unavailable norm or failed recovery is reported as a successful zero.

The test suite covers explicit sample, fit, test and overlap deficiencies; complex conjugation; physical units and coefficient normalization; observation/column permutations; test-row phase/scaling/permutation; zero fit weights; active sampled-zero columns; inactive F0/G0; rectangular/empty batches; nonmodal physical duals; and represented extended-tail rejection.

Finite input matrices can still require unrepresentable intermediate coordinates. A tested D with about `1e318` relative scale range returns `unrepresentableIntermediate` in diagnostic mode with no R and with its original B/specification preserved; strict mode raises an error. No fallback is selected. Normalized physical duals outside floating-point range are unavailable independently of sample recovery. Full matrix D storage remains quadratic in observation count; implicit/structured whitening is not implemented here. Algebraic weak recovery without D avoids observation-sized square temporaries.

No local/uncommitted DD manuscript changes, research archives, WVM consumer changes, package snapshot changes or dependency/version changes are included. Broader family adapters and physical qualification remain in #48/#58/#60; WVM adoption remains #59.
