# Neutral-label API verification

Implemented and audited against the neutral-label structure, generalized-energy endpoint conditions, and APV projection in `literature/ape-apv-free-surface/main.tex`, starting from InternalModes V2 commit `9cf40ef`. The constant-label convention also follows the manuscript's neutral-layer relabeling note.

## Scientific and API audit

- Checked analytical DD, DN, ND and positive NN profiles against the manuscript's depth normalization, including finite nonzero and infinite endpoint accelerations.
- Required explicit exclusion of the APV-invisible constant label for a fully neutral NN column. At zero Coriolis frequency, profile evaluation remains available and APV projection rejects the degenerate problem.
- Preserved the continuous APV normalization rather than fitting a sampled Gram inverse. Row scaling gives the prescribed projector an identity coefficient system while retaining the physical enstrophy Gram as its error metric. This avoids false rank loss at small nonzero rotation and does not claim grid conditioning or automatic recovery qualification.
- Checked collection column identities, repeated/chunked evaluation, selected columns, and unavailable generic Galerkin capabilities. Interval-interior second derivatives explicitly omit distributional interface terms.

## Focused MATLAB checks

All 30 distinct test methods passed across the focused runs:

| Test class | Passing methods |
| --- | ---: |
| `IMNeutralLabelModesTests` | 11 |
| `IMBasisCollectionTests` | 7 |
| `IMProjectionOnGridTests` | 5 |
| `IMBulkBasisEvaluationTests` | 7 |

The first ten-method neutral run passed nine methods and exposed an empty-chunk evaluation failure when no sample intersected an interval. Skipping empty interval rows fixed it; only the failing `intervalSupportLabelsAndCollection` method was rerun, and it passed.

The independent audit added `nonmodalAPVMatchesIndependentQuadrature` to the permanent test suite. It projects a complex polynomial APV source outside the retained span using a 64-point Gauss-Legendre rule and compares against independently integrated analytical profiles. Seven configurations cover the endpoint combinations, finite/infinite endpoint parameters, and rotation rate `1e-12`; relative coefficient errors must remain below `2e-12`. This method and the 19 existing integration regressions all passed in the final numerical batch.

All MATLAB commands used `matlab -batch` outside the sandbox, as required by the workspace. MATLAB emitted an existing startup-path warning for `/Users/jearly/Documents/MATLAB`; it did not prevent any test from running. No full-suite, coordinate-study, or production WVM qualification is claimed.

## Static and documentation checks

Code Analyzer checked all eight changed MATLAB files once. The new API and test files have no findings. The only finding is the existing `STISA` advisory on the unchanged exact-class check in `IMBasisCollection/evaluateChunks.m`; it was verified against `HEAD` and left unchanged.

Documentation links and supported inline documentation tokens passed checks. This repository has no configured `docs:check` task; the focused source checks were used instead. No website build or generated website update was performed. Whitespace and repository scope were checked; package metadata, dependencies, released snapshots, and coordinate-study archives are unchanged.
