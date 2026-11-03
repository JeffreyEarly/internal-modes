# Neutral-label API verification

The qualification below records the historical generalized-energy-dependent endpoint policy and is preserved as an implementation ledger. The current manuscript fixes both interval endpoints and excludes the constant label; the updated qualification appears first.

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

The historical grid-recovery counts below used the old Neumann surface policy for a neutral cap. A subsequent [shared-grid DD requalification](NeutralLabelGridRequalification.md) reran 38 grid/budget cases, preserved all archived counts and boundary decisions, and directly qualified 151 physical plus two DD neutral amplitudes at the exact 65-point default. Its 23 focused study tests and independent numerical audits are separate from the 70 provider test methods above. No full-suite, clean exported-package or production WVM qualification is claimed. Earlier Neumann and constant-label claims below describe only the historical API.

Implemented and audited against the neutral-label structure, generalized-energy endpoint conditions, and APV projection in `literature/ape-apv-free-surface/main.tex`, starting from InternalModes V2 commit `9cf40ef`. The original qualification below used explicit constant-label exclusion. The subsequent creation-API update retains that label automatically and gives its APV degeneracy an explicit projection error.

## Automatic-constant update

All 13 methods in `IMNeutralLabelModesTests` passed after adding automatic constant inclusion and the `"depthRMSG"` normalization name. The two additional methods check:

- A fully neutral NN column with `nModes=0` returns the constant alone, with exact unit G norm, zero derivatives, zero APV, and an explicit unobservable-amplitude error on APV projection. Zero positive-mode counts on other intervals also remain valid.
- A basis containing the constant preserves positive-mode identities and recovers selected positive amplitudes using independently formed analytical APV profiles, including the `IMBasisCollection` path. Any selected constant column is rejected rather than silently assigned zero.

Explicit `constantLabelPolicy="reject"` and `"exclude"` remain covered for compatibility. The focused run used `matlab -batch` outside the sandbox. MATLAB emitted existing package-path warnings concerning the installed InternalModes path; the tests completed successfully. This run establishes the neutral-label baseline for the creation update and does not claim reruns of the integration classes listed under the original qualification.

Whitespace checks passed on the four changed neutral-label MATLAB files. The subsequent creation-API integration removed the small local-mode-number allocation advisory; the constructor now has no Code Analyzer findings. The affected constant-label cases passed again. See [mode creation verification](ModeCreationAPIVerification.md) for the complete integration ledger.

## Scientific and API audit

- Checked analytical DD, DN, ND and positive NN profiles against the manuscript's depth normalization, including finite nonzero and infinite endpoint accelerations.
- The original implementation required explicit exclusion of the APV-invisible constant label for a fully neutral NN column. The current default retains it and requires explicit observable-column selection for APV projection. At zero Coriolis frequency, profile evaluation remains available and APV projection rejects the degenerate problem.
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
