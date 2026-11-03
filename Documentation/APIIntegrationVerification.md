# V2 API integration verification

## PR review corrections

The follow-up to `c2f7b50` resolves scalar spectral RMS products with the same oversampled physical-depth quadrature used for paired F/G profiles. Collection `projectionOnGrid` now enforces the selected basis's physical scalar capabilities; `assess` retains auxiliary EVP pairings. The neutral-label quick start supplies `variable="APV"`, and the earlier coordinate-study result below explicitly identifies its working-tree dependencies.

All **111 current cases in ten affected provider suites** are covered by passing results: 110 passed the initial focused run, and one old test that expected the collection capability bypass passed after migration to the assessment API. Five new regressions cover an exactly represented high-degree Chebyshev profile, a near-capacity scalar solve, physical-depth normalization in stretched coordinates, wave capability parity across backends, and collection page/column selection. RMS checks use an analytic integral and independent adaptive integration. The two affected class references were generated and checked in temporary output.

This follow-up does not rerun the full provider or consumer suites. The paired F/G quadrature algorithm is unchanged. High-degree RMS verification does not qualify the separate native-quadrature scalar Gram calculations. Execution evidence is in `/private/tmp/pr68-comment-fixes`.

## Post-review normalization and count hardening

Completed on 2026-11-03 from provider `6ba7b6e` and WVM `c4621b3e`. This follow-up implements the revised contracts in [ModeCreationAPI](ModeCreationAPI.md), [BasisNormalization](BasisNormalization.md), [CommonScalarProjection](CommonScalarProjection.md) and [HydrostaticRootTransform](HydrostaticRootTransform.md).

- Fixed-grid hydrostatic adapters and root-transform physical records now use the actual full normalized F/G Gram matrices and majorants. They no longer infer a metric from a replaceable name or infer the represented G metric from F/h. A finite-difference regression checks physical forms and independent off-span inputs. Root qualification, quadrature fitting and truncated-collocation recovery retain their existing policies.
- Built-in profiles remain dimensionless. The catalogue uses `rmsF`, `rmsG`, generic `rms`, `maxAbsF`, `maxAbsG`, `surfaceF`, `waveEnergy`, `generalizedPotentialEnstrophy`, `meanDensityAnomalyEnergy` and `endpointResponse` where supported. Canonical names are reserved; custom rules remain replaceable and validate divisors before MATLAB numeric conversion. Retired aliases, enum inputs, `innerProduct*` variants and `generalizedEnergy` normalization are rejected. Energy matrices and generalized-energy orthogonalization remain available.
- Bulk counts align with every requested wavenumber, including zero. Unequal repeated-wavenumber counts share one maximal solve, with exact selected prefixes and per-request provenance. Interval families require `nModesPerInterval`; physical backgrounds and applicable raw factories require explicit surface choices. Empty page and column shapes are predictable.
- `projectionOnGrid` defaults to physical preparation alone. Physical application uses `projectPhysical`/`physicalDualMatrix`; explicit recovery uses `recovery.transformForward`/`recovery.forwardMatrix` or `scalarTransform`. The ambiguous `IMProjection.project` and `IMProjection.forwardMatrix` aliases are removed. Low-level direct matrix construction retains its documented recovery default for assessment callers.
- Analytical projection provenance now reports the actual Gauss–Legendre integration rule and original-inventory point budget. Selection preserves that budget.
- WVM requests provider `waveEnergy` and applies its own F/G coordinate factors, including derivatives and k=0. Existing model coefficients, physical fields and saved-array conventions are preserved. MDA keeps its original scale under the new name. Active coordinate-study callers use explicit physical APV duals, recovery and counts. Manuscript source required no edits.

### Verification

At `c2f7b50`, the provider inventory contained **737 distinct cases**, all covered by passing results. One complete run exercised 731 cases: 727 passed and four stale test expectations failed. The corrected four passed focused reruns. Six additional regression cases passed separately: five custom-divisor cases and the finite-difference full-Gram case. Five targeted checks also passed after the final root/provenance corrections. The final inventory comparison found no unverified cases. This is a full run plus targeted repairs and additions, not a second full-suite run.

| Check | Result |
| --- | --- |
| Provider V2 | 737 distinct cases covered; no unresolved failures |
| Affected WVM suites and field/source/restart checks | 30 passed |
| Affected coordinate-grid study suites | 31 passed in the working tree, including pre-existing uncommitted DD neutral-label fixes; not independent qualification of clean companion commit `f459aa8` |
| Independent neutral-label validation | 154 checks passed |
| Original-code NetCDF restart | Exact arrays, coefficients, fields and energy |
| Two-profile old/new WVM comparison | 66 quantities; maximum relative difference 4.25e-15 |
| Old/new F/G and spectral derivatives, positive k and k=0 | Maximum relative difference 4.395e-15 |
| Manuscript source-array smoke | All six scripts passed at original resolutions |
| Bounded study/tool smoke | 9 reference families, 38 prefix rows; one constant-profile calibration |
| Code Analyzer | No new findings after correction; existing advisories retained |
| Scope, whitespace, manifest and Markdown links | Passed; no generated website or package payload changes |

Provider execution evidence is in `/private/tmp/im-v2-implementation`; consumer evidence is in `/private/tmp/im-v2-caller-migration-20261103`. Temporary runner-path errors were corrected without changing source behavior. No `docs:check` target is configured in this provider checkout; canonical Markdown links and source documentation were checked, and website generation was excluded by the requested scope. No local asset prevented completion.

Existing user edits in the provider example and coordinate study were preserved. The coordinate-study tests depended on those existing DD fixes, including removal of `constantLabelPolicy` and updated endpoint expectations. Those fixes remain outside companion commit `f459aa8`; the 31 passing cases describe the tested working tree, not that commit alone. Full performance campaigns and archived grid-qualification datasets were not regenerated. A complete non-surface-active `WVTransformBoussinesq` port, model-owned grid calculus choices and release/export remain separate work. This verification preceded the subsequent commits and provider push.

## Earlier integration record

This record describes the earlier integration. For the current normalization, exact-count and projection contracts, see [ModeCreationAPI](ModeCreationAPI.md), [BasisNormalization](BasisNormalization.md) and [CommonScalarProjection](CommonScalarProjection.md). Historical test counts and retired API names below are preserved as evidence of that earlier state.

Implemented and verified on 2026-11-03 for PR [#68](https://github.com/JeffreyEarly/internal-modes/pull/68). The integration merges `InternalModesEVP` at `89ab5a7501189d1bf11adb471a03b03e9781719a` into the capable-basis migration at `f144924d72b1250f053e1886201e7c4880616926`. All 23 assessed source/test/documentation conflicts are resolved. See the [implementation plan](APISimplificationPlan.md) for contracts and scope.

## Delivered behavior

- Construction returns capable family bases through `IMBackground`, family problems and explicit solver requests. The removed namespace, wrapper and old background class stay removed.
- Structured capabilities, cheap scientific/construction records, horizontal identities and lazy normalization are shared across numerical, analytical and specialized bases. Selection preserves original construction evidence and integration effort while validating only retained normalization factors.
- Counts remain exact totals. Named evaluation, per-call normalization, explicit-empty columns and prevalidation before collection callbacks remain consistent.
- Boundary-buoyancy modes use explicit endpoint inventories, `orthogonalizeGeneralizedEnergy(g0=...,gd=...)`, signed `boundaryDepth` and `boundaryResponseMatrix`. General public rotations are removed. Orthogonalization requires the complete originally chosen endpoint family; empty and one-endpoint constructions remain valid. Negative, zero and repeated boundary depths are retained.
- Endpoint scientific normalization names `endpointResponse` and `generalizedEnergy` cannot be replaced by custom registration. Custom and profile-rule names remain replaceable. This prevents silently changing canonical geometry or null-energy meaning.
- Native scalar physical adapters retain independent source observations and physical duals separately from corrected recovery. Generic matrix recovery, strict rank qualification and failed-recovery diagnostics remain available.
- Hydrostatic root transforms reuse captured full-span samples when binding selected prepared recovery, avoiding additional factorization and BLAS-shape roundoff mismatches. Galerkin and truncated-collocation policies remain distinct.
- Prescribed sampling and its dedicated class, example, tests, fixture, documentation and orphan support are removed. Ordinary basis evaluation supplies point observations.

## Provider tests

The final V2 inventory contains **723 distinct test cases**, all covered by passing results. The complete run exercised 721 cases and initially passed 704. Its 17 failures were remaining test callers that omitted the required projection variable or expected a retired validation error identifier. Those callers were migrated and all 17 failing cases passed their focused rerun. Two additional endpoint custom-rule regressions passed after final review. The final test-inventory comparison confirms no missing cases, unresolved failures or incomplete results. This is one full run followed by focused repairs and additions, rather than a second full-suite run.

Earlier focused integration gates covered 112 native/evaluation/count cases, 88 recovery/projection/collection/root cases, 57 endpoint/composition cases and 44 capability/identity/provenance/normalization cases. Their counts overlap with the full inventory and are not added to it.

Scientific tolerances were preserved. Focused verification caught and repaired extrema normalization, analytical phase orientation, selected custom-vector indexing, normalization cache detachment/reuse, retained-only admissibility, built-in/custom replacement precedence and nonreal factors. Endpoint review additionally repaired lazy null-energy validation, selected custom-rule validation and protection of the two scientific normalization names.

Endpoint checks include canonical unit buoyancy responses, coupled-equation residuals, coefficient orthogonality, signed and null generalized energy, repeated eigenspaces, current-coordinate form congruence, selection restrictions, absent-endpoint infinite coefficients, complex independent sources and physical/recovery covariance. Recovery checks distinguish inactive hydrostatic G0 from active constant MDA G, preserve supplied points/weights and test both root policies.

The migrated neutral-label verification tool and all four affected endpoint/composition examples pass at their defaults. The analytical exponential endpoint example agrees with its refined numerical reference within 2.3e-15 relative F error and 7.3e-16 relative G error.

## Manuscript figures

All seven mode PDFs were regenerated in a temporary copy of the manuscript MATLAB tree. All **28 saved plotted arrays are exactly equal** to their preserved pre-integration values. All seven regenerated pages are **pixel-identical at 150 dpi**, and the rendered pages were visually inspected. There are no scientific or visible figure changes.

The comparison covers Figures 2–6 and the two constant-stratification geostrophic plots. The coordinate-transformation illustration is independent of the mode API and was not regenerated. Existing manuscript sources and committed PDFs need no further edits; the earlier migration remains on manuscript `master` at `82fb249`.

## WaveVortexModel

All **23 bounded consumer test cases pass**: scalar/bulk construction, physical reconstruction, energy and source pairings, forced evolution, restart, stored-mode restoration, reuse, convergence preparation, per-wavenumber assessment and independent shooting.

Fresh constant- and exponential-stratification snapshots compare 66 numeric quantities. Profiles, spectra, reconstructed fields and energies are unchanged; three source-projection quantities differ only by rounding, with maximum relative L2 error **1.303e-15**. The provider wave relation remains F=D Gz, and the existing consumer h/D conversion is applied exactly once.

No further consumer source edits were necessary. The existing WVM migration remains on `feature/im-v2-api-migration` at `c4621b3e258be159fc9cff101203118a53a7f764`. Its repository artifact policy check passes. The complete WVM suite was not run.

## Static and documentation checks

Code Analyzer covered all **110 distinct changed MATLAB files**, with targeted rechecks after repairs. There are no parse/error diagnostics. Fifteen advisories were reviewed: allocation, scientific property names, exact-class dispatch, unused values/inputs and deliberate normalization-registry lookup in setters. Final repaired test and endpoint files report no additional findings.

Canonical inline API generation produced **438 Markdown pages for 15 classes** in temporary output. Required count/evaluation signatures, endpoint operation/properties and absence of retired rotation pages were checked. The two endpoint classes were regenerated after their final help edits. Markdown fences/local links, whitespace, removed callers, conflict markers, repository scope, unchanged manifests and generated-artifact scope pass. The provider has no configured `docs:check` task; temporary extraction and bounded source checks were used.

No package metadata, released OceanKit snapshots or generated website output changed. The user's uncommitted free-surface example edit was preserved separately from the integration and remains unmodified.

## Qualification limits and artifacts

This verifies authoring checkouts and the bounded integrations above. Numerical normalization quadrature and the available construction records do not certify convergence for arbitrary profiles or parameters. Provider release/export, WVM dependency-floor and CI-pin updates, clean exported-package qualification and a broader `WVTransformBoussinesq` port remain separate work. No required local asset was missing.

Execution artifacts are outside the repositories in `/private/tmp/im-v2-integration-20261103`, with endpoint-focused artifacts in `/private/tmp/im-v2-endpoint-work`. The combined inventory report records all 723 passing cases. No research outputs were added to WVM or the manuscript.
