# V2 API consolidation verification

This record describes the earlier integration. For the current normalization, exact-count and projection contracts, see [ModeCreationAPI](ModeCreationAPI.md), [BasisNormalization](BasisNormalization.md) and [CommonScalarProjection](CommonScalarProjection.md). Historical test counts and retired API names below are preserved as evidence of that earlier state.

Implemented locally on 2026-11-02 against the authoring checkout, preserving the preexisting uncommitted neutral-label changes and the preceding manuscript migration. No source push, tracker communication, package release or website publication was performed.

## Delivered contracts

- Removed `IMModes` and `IMModeBasis`; construction returns the capable numerical, analytical or specialized family basis.
- `solveEVP(...,nModes=N)` returns exactly N total columns. Explicit branches require matching count vectors; zero requests are empty and unavailable inventories fail explicitly.
- Added authoritative background-first family problems, a real analytical solver and exact repeated-page analytical/numerical bulk construction.
- Preserved lazy shared normalization, per-call overrides, custom rules, aligned selection, native diagnostics, signed pairings, projections, endpoint rotations and hydrostatic root transforms.
- Canonical zero-APV responses use explicit `endpoints`, independently of generalized-energy coefficients.
- All public wave paths use F = D Gz. Hydrostatic and balanced relations remain distinct. WaveVortexModel explicitly requests its existing G normalization and converts provider F by h/D at the consumer boundary.

See [API usage and migration examples](ModeCreationAPI.md) and the [implementation plan](APISimplificationPlan.md).

## Provider verification

The full `UnitTestsV2` run exercised 589 methods. It initially passed 557 and identified 32 failures. Targeted repairs and reruns resolved every failure, with four additional regression methods bringing the combined inventory to **593 distinct passing methods**, zero unresolved failures and zero incomplete results. This is a full run followed by focused reruns, not a second full-suite run.

Scientific tolerances were retained. Most assertion changes reflect the agreed count, normalization-name, metadata and derivative-capability contracts. The full run also exposed an inherited smoothed-bandwidth test assumption: classical N and rotation-aware bandwidth have distinct coordinate scales even when both produce the same uniform physical grid. The archived baseline reproduced that failure; the test now checks each independent defining formula without relaxing its bounds.

Integration review found and repaired hydrostatic auxiliary solves that omitted the retained barotropic column. Additional regressions cover root generation and leakage extension, custom diagnostic derivative rejection, stored companion derivative consistency, and custom normalization/dual preservation under selection. Arbitrary diagnostic callbacks no longer receive an assumed hydrostatic derivative.

## Manuscript figures

All eight supporting MATLAB scripts use the canonical API. The seven manuscript PDFs were regenerated and visually inspected. Every rendered pixel at 150 dpi matches the preserved pre-consolidation PDF baseline; there are **no visible figure changes**.

Across 28 saved plotted arrays, all frequencies, equivalent depths and eigenvalues are bitwise unchanged. The largest normalized-profile difference is 4.44e-16; the largest penetration-depth difference is 4.55e-13 m. The nonlocal endpoint experiment remains mathematically unchanged. The direct frequency-squared study retains its independent pencil; its frequency and depth outputs agree with its previous recorded results within 4.2e-15 relative error.

## WaveVortexModel verification

Twenty distinct bounded consumer methods passed: bulk/scalar construction, optimized convergence preparation, physical reconstruction, energy, source pairings, forced evolution, restart, stored-mode restoration, reuse checks, per-wavenumber assessment and independent shooting qualification. The previously recorded shooting limitation did not reproduce: the unchanged focused method passed with both archived and current provider/consumer sources.

Separate fresh MATLAB processes compared archived and current physical states for constant and exponential stratification. All 66 compared numeric quantities agree within 1.03e-15 relative L2 error. Wave G, its physical forward operator, depths, frequencies and inertial depths are bitwise equal; velocity, pressure, energy and source projections agree to rounding.

## Static checks and scope

Code Analyzer covered 92 changed MATLAB files; the two files receiving final naming/help repairs were rechecked. Remaining findings are advisory naming, allocation, exact-class, unused-input or existing exploratory-code diagnostics, with no parse or error findings. Canonical inline API documentation was generated into a temporary directory for 15 classes (355 Markdown pages). Markdown fences and local source links, whitespace, runtime references to the removed classes, package manifests and repository artifact scope were checked.

The provider has no configured `docs:check` task. Generated website content, package versions/dependencies and released OceanKit package snapshots are unchanged. The WaveVortexModel artifact policy check passes. Verification outputs and baseline archives are outside the checkouts in `/private/tmp/im-api-consolidation`; no recorded study outputs were added to WaveVortexModel.

This handoff qualifies the authoring checkouts and the bounded integrations above. A clean exported-package release and the complete WaveVortexModel suite were not run. No required local data asset was missing.

## Background naming cleanup — 2026-11-03

Renamed the physical-column class to `IMBackground`, including its constructor and source factories, type checks, error identifiers, callers, tests and documentation. The mechanical rename was checked against preserved pre-edit source in 36 files across the provider, WaveVortexModel and manuscript. No scientific equation, default or figure output changed.

The 50 focused provider methods and two consumer construction/persistence methods passed. Code Analyzer covered the 30 renamed MATLAB sources, with only two existing unused-input advisories in the analytical solver. Inline API documentation for nine affected classes was generated into a temporary directory. Whitespace, documentation fences/links, removed runtime references, unchanged manifests/generated-site scope and the consumer artifact policy passed. The broader numerical and figure comparisons above remain applicable; figures were not regenerated for this naming-only change. Rename verification artifacts are in `/private/tmp/im-background-rename`.

## API usage consistency — 2026-11-03

Numerical and analytical `solveEVP` calls now require an explicit `nModes`; positive bulk pages do too, while zero-only requests use `nInertialModes`. Scalar, specialized and collection evaluators use explicit `variable` and optional `derivativeOrder` name-value arguments. Ordinary and streamed collections accept per-call normalization without changing stored defaults. Collection capability metadata reflects each basis's supported derivatives, and unsupported requests fail before streaming callbacks. Endpoint normalization remains tied to an explicit boundary-coordinate rotation.

All **127 distinct focused test methods pass**, including fourteen new count/evaluation methods and the affected existing test classes. Tests cover numerical, finite-difference and analytical backends; empty inventories; custom normalization; native and fallback evaluation; repeated pages and columns; chunking; nonlinear companion recovery; and endpoint conventions. A finite-difference interpolation bug exposed by the zero-count tests was repaired so empty profiles and derivatives retain sample-by-zero shapes.

Two older test assumptions were corrected. Generic scalar bases now support `evaluate`, so their test checks the explicit `u` evaluation contract. The generalized-enstrophy derivative test's absolute algebraic-equality assertion failed identically against the preserved baseline. These companion derivatives are differentiated independently from native coefficients; the test now uses the established geostrophic per-column differential-equation residual criterion of `1e-7`. Its measured maximum dimensionless residual is `7.4e-10`. No derivative algorithm was changed.

All four affected examples pass at their full defaults: free-surface exponential waves, unified basis assessment, wave counts by wavenumber, and the bandwidth-coordinate comparison. The corrected free-surface boundary residual is `7.0e-16`. The bandwidth example now explicitly selects `"bandwidth"` for the sqrt(B) coordinate. All twelve nonneutral archived comparison counts agree; two neutral results now report an unavailable exact inventory with `NaN` counts and the solver diagnostic, preserving the requests of 62 and 126 columns rather than silently returning 59 and 119. Manuscript figure sources and PDFs were not changed or regenerated in this batch.

Code Analyzer covered all 38 changed MATLAB files. Its 33 advisories occur on unchanged baseline source lines; no parse errors or new findings were introduced. Inline API documentation was generated into a temporary directory for the changed public classes. The repository has no configured `docs:check` task; bounded generated-content checks and canonical Markdown fence/link checks were used. Package metadata, dependencies, released snapshots and generated website files are unchanged. Verification artifacts and the pre-edit source baseline are in `/private/tmp/im-api-consistency`; no required local asset was missing.
