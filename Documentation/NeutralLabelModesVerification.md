# Neutral-label API verification

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
