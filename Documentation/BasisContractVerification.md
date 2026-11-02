# Common basis contract verification

Verified on 2026-11-02 from `InternalModesEVP` at `16235ef` (the merged manuscript DD correction). This increment implements #44–46 without changing the recovery engine, observation-grid design or WVM code.

All MATLAB commands used `matlab -batch` outside the local sandbox, with one computational thread. Runtime: MATLAB R2026b Update 1, MACA64.

## Independent review

The provenance/factory and capability/collection implementations had separate authors. Each reviewed the other implementation; a third agent performed read-only mathematical, numerical and framework review. Root-owned constructor, horizontal identities and column selection were also independently reviewed.

The independent runtime audit passed **310 controls across all 11 families**. Controls included:

- Fixed-frequency signed dispersion above/below the inertial frequency, both rotation hemispheres, and negative equivalent depths above the buoyancy frequency.
- Scalar/collection evaluation agreement for declared derivatives and aliases, repeated/chunked pages and empty/reordered inventories.
- Scientific versus construction records, active column identity, physical/generalized-energy boundary separation and DD semantics.
- External APV page wavenumbers, genuine inertial k=0 and rejection of contradictory horizontal page assignments.
- Requested-empty generalized-enstrophy default normalization and parity with selected-empty capabilities.

The review found and corrected an empty-basis normalization regression: inspecting only the retained source payload lost generalized-enstrophy capability when the backend returned no columns. The shared capability builder now reads the existing canonical EVP metadata without solving, quadrature or projector preparation. Explicit intrinsic F-prime declarations also preserve the derivative contract for analytical barotropic and empty endpoint sectors. Independent requalification passed. No substantive mathematical, numerical or framework finding remains open.

## Unit tests

| Coverage | Distinct focused methods |
| --- | ---: |
| New capability, scalar/collection parity and empty-inventory controls | 10 passing |
| New requested/resolved scientific and construction provenance | 9 passing |
| Existing common and lower-level collection controls | 15 passing |
| New horizontal identity and exact rigid-lid dispersion controls | 6 passing |

Related source changes were batched; failing and affected author controls were rerun selectively. Three initial horizontal-test failures came from test fixtures/assertions (unsupported analytical MDA, expected spectrum ordering and 0-by-0 versus 1-by-0 empty shape); their corrections are covered by the final integration run.

Full V2 integration at `61981d7`: **645 of 647 methods passed**. All 25 new methods passed, including all six horizontal-identity controls; all 40 distinct affected focused methods passed. Production source was frozen before this run and before the consolidated Code Analyzer pass. The subsequent test-fixture cleanup below resolves the two failures.

The predecessor's full run passed 620 of 622 methods. The `61981d7` run had the same two known baseline failures: `IMEigenvalueProblemRefactorTests/developerBasisMethodsArePublicAndCallable` (analytical/EVP stratification mismatch; this method was also incomplete) and `IMSmoothedBandwidthCoordinateTests/constantControlKeepsTheGridAndPhysicalSpectrum` (an incorrect classical-WKB monitor comparison). Their implementations and fixtures were unchanged at that revision. No new full-suite failure was introduced.

## Baseline test-fixture cleanup

The developer-method test now gives the analytical solution a matching constant-stratification EVP, retaining the varying-stratification numerical fixture and all existing method checks. The rotating constant-monitor test now compares smoothed bandwidth with unsmoothed bandwidth, and independently checks the exact monitor at latitude 45 degrees. Grid, physical-stratification and eigenvalue checks are preserved, with unchanged tolerances. Production code is unchanged.

Independent mathematical/framework review found no actionable issue. Both corrected methods passed, together with the rotating filter-residual control and the constant-solution tests that reject different/interior-mismatched stratification: **5 of 5 methods passed**, none incomplete. These passing results resolve both failures in the earlier full run. The full suite was not repeated for this test-only cleanup.

Cleanup Code Analyzer checks found no new messages: the smoothed-bandwidth test file has none, and all 14 advisories in the developer-method test file match the preceding revision. The latter file was compared with its baseline to establish that the reported messages were pre-existing. Source documentation and scope checks passed; test-only changes require no API-documentation regeneration.

## Static, documentation and scope gates

- Code Analyzer ran once on all 13 changed MATLAB files. No new findings; one existing exact-class advisory remains on an unchanged native-grouping statement in `IMBasisCollection/evaluateChunks.m`.
- Temporary API documentation extraction succeeded once for `IMModes`, `IMModeBasis` and `IMBasisCollection`, producing 56 Markdown files. The new public record and horizontal-identity entries were present.
- Source Markdown local links and affected fences were checked. Existing legacy README fences were preserved. There is no configured `docs:check` task in this repository; source checks and temporary extraction supply the relevant documentation gates.
- Whitespace, repository scope, unchanged package manifest/dependencies and unchanged generated website artifacts were checked. No released OceanKit snapshot was edited. The unrelated original dirty manuscript worktree was preserved byte-for-byte in its tracked diff.

## Scope of the evidence

Construction counts and backend settings remain explicitly unverified solve/convergence evidence. Capability records describe available operations and conditional preparation requirements; they do not certify rank, physical completeness, noise behavior or an inverse on a particular grid. Scientific pairing availability does not imply an implemented normalization operation: the analytical barotropic-only bypass has G=0 and cannot use generalized-energy G normalization. Common pairing-based normalization remains a separate increment (#47). This increment adds no new physical projection adapter, normalization-changing operation or WVM adoption. It preserves externally imposed observation grids as a supported input to the shared recovery layer.
