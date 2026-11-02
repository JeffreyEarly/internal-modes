# Neutral-label API verification

The qualification below records the historical generalized-energy-dependent endpoint policy and is preserved as an implementation ledger. The current manuscript fixes both interval endpoints and excludes the constant label; the updated qualification appears first.

## Reviewed delivery on the shared recovery stack

Integrated on 2026-11-02 against merged shared-recovery revision `91cc20ce3cb7fe386639c20307a753e62e2264fa`. The reference is `literature/ape-apv-free-surface/main.tex`, SHA-256 `eef2e7f2ff899557d6eaba0030116df7f58814f6f61892b3ef486ea1df41958c`: equations `eq:neutral-label-vertical-structure`, `eq:neutral-label-simple-boundary-conditions`, `eq:neutral-label-apv-normalization`, `eq:projection-neutral-label-functionals`, and Appendix `appendix:neutral-label-projection`. The governing formulas and migration are reproduced in the [API guide](NeutralLabelModes.md); verification does not require a private manuscript or research cache.

The initial port passed all 70 focused methods in the seven classes listed below. Independent mathematical, numerical and framework reviews checked code they did not author. A framework review found an empty selected common view preparing discarded APV norms before taking a zero-column prefix. It now constructs only an empty analytical source and uses the normal grid validation path. This preserves zero-column projection at tiny nonzero rotation even when the discarded modes' APV norms underflow. The new regression, the existing selected/empty regression and a new physical-versus-recovery regression all passed; the framework reviewer independently reran both new methods and an additional shifted, reordered `maxAbsG` empty view.

The prescribed neutral APV dual remains a physical `IMProjection`, with `projectPhysical == project`, no recovery object, and unsupported Gram/round-trip assessment. A three-node control deliberately gives physical `B*A = diag([1 2])`, while separately requested shared-engine sample recovery gives `R*A = I`. Quadrature error is retained rather than silently corrected. Continuous positive APV norms are not sampled-grid observability claims.

The portable [validateNeutralLabels.m](../tools/validateNeutralLabels.m) publishes the numerical auditor's 154 controls, using a shifted physical domain, separated bottom/interior/surface intervals and a full-column case. Independent Gaussian/adaptive references test sine derivatives, APV and depth norms, complex nonmodal sources, both source endpoint terms, affine cancellation, both Coriolis signs, rotation rates down to `1e-100`, normalized/reordered collections, empty shapes and exact zero-depth-wave energy. An unresolved-rule control keeps its physical projection error. The largest nonzero-reference relative error was `5.89e-14` for the displacement-source integration-by-parts reference; zero-reference metrics are absolute. A separate mathematical audit obtained APV Gram relative error `2.57e-15` and source-functional absolute difference `3.56e-13` on another shifted domain.

Run with declared package dependencies on the path:

```matlab
maxNumCompThreads(1);
addpath('UnitTestsV2','UnitTestsV2/TestSupport','tools');
report = validateNeutralLabels;
results = runtests({'IMNeutralLabelModesTests','IMModeCreationTests','IMModeBasisTests','IMModeBasisCollectionTests','IMProjectionOnGridTests','IMBasisCollectionTests','IMBulkBasisEvaluationTests'});
```

The portable tool restores caller path/RNG and leaves thread configuration unchanged. The implementation changes no grid construction, quadrature weights, retained-count default, shared recovery policy, dependency/version or package snapshot. Local grid-requalification and transform-investigation reports remain outside this bounded delivery; current cross-grid retained-count claims require #58's portable qualification.

Current verification used MATLAB `26.2.0.3433840 (R2026b) Update 1`, `MACA64`, one numerical thread, with the declared package dependencies from authoring checkouts. The initial 70 methods plus the two new methods give 72 distinct passing focused methods. The portable 154-control tool also passed after its fixture background was refined to keep positive stratification between the declared maximal neutral intervals. Code Analyzer reported zero findings across all 11 changed MATLAB files; only the tool was rechecked after that refinement. Temporary API generation produced 74 Markdown pages for the four affected classes. Affected local links, supported inline tokens, whitespace, scope, unchanged package/fixture manifests and unchanged generated artifacts were checked; existing README fence formatting outside the edited lines was preserved. There is no configured `docs:check` task. No released-package or production WVM qualification is claimed.

The final full V2 run passed **620 of 622 methods**. Both failures are the same baseline defects previously reproduced on the exact clean predecessor during shared-engine delivery: `IMEigenvalueProblemRefactorTests/developerBasisMethodsArePublicAndCallable` errors on an analytical/EVP stratification mismatch, and `IMSmoothedBandwidthCoordinateTests/constantControlKeepsTheGridAndPhysicalSpectrum` fails its smoothed-bandwidth reference comparison. Their implementations, fixtures and assertions are unchanged here. Every neutral-label and shared-recovery test passed; no new full-suite failure was introduced.

## Fixed-endpoint manuscript alignment

Verified locally on 2026-11-02 against the manuscript's neutral-label structure and projection equations. Every interval uses DD sine modes, including at physical boundaries; zero, finite signed and infinite generalized-energy coefficients do not alter profiles. There is no constant-label sector. The direct constructor retains `g0`/`gd` as provenance, and the common factory accepts optional `energyBoundary` provenance. The former `constantLabelPolicy` option is removed.

All 70 distinct methods passed across the following focused suites:

| Test class | Passing methods |
| --- | ---: |
| `IMNeutralLabelModesTests` | 15 |
| `IMModeCreationTests` | 17 |
| `IMModeBasisTests` | 12 |
| `IMModeBasisCollectionTests` | 7 |
| `IMProjectionOnGridTests` | 5 |
| `IMBasisCollectionTests` | 7 |
| `IMBulkBasisEvaluationTests` | 7 |

The first batch passed 67 of 68 methods. The new affine-source cancellation assertion differed from zero by `1.67e-16` with a `1e-16` tolerance. Its tolerance was changed to `5e-15` to allow floating-point summation error; only that failing method and two subsequently extended creation methods were rerun, and all passed. No numerical source correction was needed for that roundoff-level test failure.

Independent formulas in the tests check sine profiles and derivatives for interior, surface, bottom and whole-column intervals; exact endpoint zeros; one-sided derivative traces and zero extension; depth and enstrophy normalization; zero counts and scalar count broadcasting; complex and nonmodal APV recovery; interval-constant APV orthogonality; and the integrated-by-parts displacement-source functional with nonzero endpoint values and affine cancellation. Common-factory tests reject even `N2=1e-30` as a declared exact neutral interval, preserve scientific identities, and retain the exact zero-depth wave regression.

An independent audit then found two empty-inventory integration failures: a raw zero-count neutral basis was rejected by collection metadata validation, and an empty selected common view forwarded an empty column selection as all source columns. Both are fixed. Two new permanent regression methods cover raw and common empty inventories, repeated-page evaluation and projection, and reordered `maxAbsG` metrics. Seven new or affected methods passed after the fixes. Prescribed-dual and recipe provenance now report the selected, normalized continuous APV Gram, including a `0x0` matrix for an empty view. A final metadata check extended this guarantee to raw reordered-column projectors and their collection path. The affected selected-column regression passed after updating only provenance; the numerical dual was unchanged, and Code Analyzer remained clean for the projector and test.

The independent final audit passed all seven groups after these fixes. It used a shifted physical domain of depth 1321 m; interior, surface, bottom and full-column intervals through 13 modes; signed/infinite energy coefficients; both Coriolis signs and rotation as small as `1e-12`; independent Gauss/adaptive Grams and zero APV integrals; complex sources outside the retained span with both volume and endpoint terms; affine cancellation; disjoint, reordered and scaled common/collection projections; raw, factory and selected empty paths; and analytical critical `h=0` wave-energy formulas. The largest normalized residual was `1.14e-13` for affine cancellation; all other reported residuals were at most `6.6e-15`. These are independent audit groups, separate from the 70 permanent test methods.

Code Analyzer reported no findings across all ten changed MATLAB files. The four MATLAB files changed by the audit fixes were checked after those fixes; the other six retain their prior successful result. Whitespace, local Markdown links and fences, supported inline documentation tokens, repository scope, unchanged package manifest/dependencies and unchanged generated artifacts were checked. This repository has no configured `docs:check` task; static documentation checks were used. No website build or generated website update was performed. All MATLAB commands used `matlab -batch` outside the sandbox; existing package path-order warnings did not prevent execution.

The historical grid-recovery counts below used the old Neumann surface policy for a neutral cap. Subsequent local DD study requalification and transform-investigation reports are not versioned in this provider change. Portable cross-grid joint-recovery evidence with the current shared engine remains [#58](https://github.com/JeffreyEarly/internal-modes/issues/58). Earlier Neumann and constant-label claims below describe only the historical API, not current DD retained-count or aliasing guarantees.

The historical implementation was audited against the then-current neutral-label structure, generalized-energy endpoint conditions, and APV projection in `literature/ape-apv-free-surface/main.tex`, starting from InternalModes V2 commit `9cf40ef`. The original qualification below used explicit constant-label exclusion. The subsequent historical creation-API update retained that label automatically and gave its APV degeneracy an explicit projection error; both policies are superseded by DD.

## Historical automatic-constant update

All 13 methods in `IMNeutralLabelModesTests` passed after adding automatic constant inclusion and the `"depthRMSG"` normalization name. The two additional methods check:

- A fully neutral NN column with `nModes=0` returns the constant alone, with exact unit G norm, zero derivatives, zero APV, and an explicit unobservable-amplitude error on APV projection. Zero positive-mode counts on other intervals also remain valid.
- A basis containing the constant preserves positive-mode identities and recovers selected positive amplitudes using independently formed analytical APV profiles, including the `IMBasisCollection` path. Any selected constant column is rejected rather than silently assigned zero.

Explicit `constantLabelPolicy="reject"` and `"exclude"` were covered for compatibility in that historical API. The focused run used `matlab -batch` outside the sandbox. MATLAB emitted existing package-path warnings concerning the installed InternalModes path; the tests completed successfully. This run established the neutral-label baseline for that creation update and did not claim reruns of the integration classes listed under the original qualification.

Whitespace checks passed on the four changed neutral-label MATLAB files. The subsequent creation-API integration removed the small local-mode-number allocation advisory; the constructor now has no Code Analyzer findings. The affected constant-label cases passed again. See [mode creation verification](ModeCreationAPIVerification.md) for the complete integration ledger.

## Historical scientific and API audit

- Checked analytical DD, DN, ND and positive NN profiles against the manuscript's depth normalization, including finite nonzero and infinite endpoint accelerations.
- The original implementation required explicit exclusion of the APV-invisible constant label for a fully neutral NN column. The then-current default retained it and required explicit observable-column selection for APV projection. At zero Coriolis frequency, profile evaluation remained available and APV projection rejected the degenerate problem.
- Preserved the continuous APV normalization rather than fitting a sampled Gram inverse. Row scaling gives the prescribed projector an identity coefficient system while retaining the physical enstrophy Gram as its error metric. This avoids false rank loss at small nonzero rotation and does not claim grid conditioning or automatic recovery qualification.
- Checked collection column identities, repeated/chunked evaluation, selected columns, and unavailable generic Galerkin capabilities. Interval-interior second derivatives explicitly omit distributional interface terms.

## Original focused MATLAB qualification

All 30 distinct test methods passed across the original focused runs, before the two automatic-constant tests were added:

| Test class | Passing methods |
| --- | ---: |
| `IMNeutralLabelModesTests` | 11 |
| `IMBasisCollectionTests` | 7 |
| `IMProjectionOnGridTests` | 5 |
| `IMBulkBasisEvaluationTests` | 7 |

The first ten-method neutral run passed nine methods and exposed an empty-chunk evaluation failure when no sample intersected an interval. Skipping empty interval rows fixed it; only the failing `intervalSupportLabelsAndCollection` method was rerun, and it passed.

The independent audit added `nonmodalAPVMatchesIndependentQuadrature` to the permanent test suite. It projects a complex polynomial APV source outside the retained span using a 64-point Gauss-Legendre rule and compares against independently integrated analytical profiles. Seven configurations cover the endpoint combinations, finite/infinite endpoint parameters, and rotation rate `1e-12`; relative coefficient errors must remain below `2e-12`. This method and the 19 existing integration regressions all passed in the final numerical batch.

All MATLAB commands used `matlab -batch` outside the sandbox, as required by the workspace. MATLAB emitted an existing startup-path warning for `/Users/jearly/Documents/MATLAB`; it did not prevent any test from running. No full-suite, coordinate-study, or production WVM qualification is claimed.

## Original static and documentation checks

Code Analyzer checked all eight changed MATLAB files once. The new API and test files have no findings. The only finding is the existing `STISA` advisory on the unchanged exact-class check in `IMBasisCollection/evaluateChunks.m`; it was verified against `HEAD` and left unchanged.

Documentation links and supported inline documentation tokens passed checks. This repository has no configured `docs:check` task; the focused source checks were used instead. No website build or generated website update was performed. Whitespace and repository scope were checked; package metadata, dependencies, released snapshots, and coordinate-study archives are unchanged.
