# Grids, physical projections and discrete mode recovery

Investigation measured on 2026-11-02 against provider revision `8717aca`, including the uncommitted manuscript-aligned DD neutral-label changes. This is an investigation and implementation recommendation. The subsequent API consolidation removed the creation facade and wrapper; the responsibilities below are stated for the current capable family bases. The broader recovery-engine proposal remains separate.

The subsequent implementation plan is tracked in [shared recovery engine #57](https://github.com/JeffreyEarly/internal-modes/issues/57) and [cross-grid qualification #58](https://github.com/JeffreyEarly/internal-modes/issues/58). The user approved a clean v2 API redesign and provider completion before production WVM adoption, tracked separately in [follow-up #59](https://github.com/JeffreyEarly/internal-modes/issues/59). Compatibility suggestions below are historical recommendations rather than implementation constraints: preserve scientific meaning, and deliberately replace obsolete scalar-only constructors/diagnostic fields as the new contracts require. The issue plan also distinguishes a least-squares fitting factor from the positive observation/noise norm and makes the active-subspace algebra explicit.

The stack can have one coherent architecture across root, depth, mapped, smoothed-bandwidth and custom grids. The essential distinction is between continuous modes, observations of those modes, physical projection functionals, and recovery from finite samples. A grid selects sampling geometry and quadrature. It does not uniquely specify a forward transform. The current hydrostatic root implementation already demonstrates this distinction, while the common fixed-grid bridge and assessment APIs provide useful foundations. General full-state recovery and compatible grid calculus are still incomplete.

## The objects and their responsibilities

| Layer | Current owner | Contract to preserve or complete |
| --- | --- | --- |
| Physical problem and EVP | `IMBackground`, `IMGeneralizedEnergyBoundary`, `IMBoundaryCondition`, `IMInternalModes` | Stratification, domain, rotation, physical surface convention, independently specified generalized-energy endpoint coefficients and the differential eigenproblem. |
| Continuous solution | `IMSolver`, capable numerical/analytical and specialized family bases | Scientific column identities, normalization, continuous evaluation, supported derivatives/traces, physical coefficient forms and projection recipes. Numerical EVP resolution is independent of physical sample count. |
| Sampling geometry | `IMVerticalGrid` | Ordered physical points, quadrature weights, domain and design provenance. Root points are a grid. Solver-native points are also a possible grid. |
| Observations and calculus | Currently distributed between provider recipes and WVM/study code | Define the sampled variables, independent boundary observations, reconstruction/interpolation, derivatives/integrals/traces and data units. Declare the positive observation norm separately from physical quadrature. |
| Prepared transform | Scalar `IMDiscreteTransform`; aligned F/G `IMInternalModesDiscreteTransform`; study joint matrices | Sampled synthesis, separately retained physical dual, explicitly chosen recovery policy and its diagnostics. General complex, coupled observations need one reusable numerical engine. |
| Qualification and selection | Assessment measurements plus caller decisions | Independently check solve accuracy, quadrature, recoverability, physical source projection, noise, unresolved inputs, boundary shapes and products. Keep measurements separate from acceptance thresholds and band selection. |

Family problems and solvers describe continuous-family creation. A named combined convenience factory can solve a basis, design a grid and prepare a transform, provided it returns the reusable continuous basis and records each stage. `IMInternalModesDiscreteTransform.hydrostaticRoot` does this. It should share the fixed-grid recovery engine rather than establish a parallel hierarchy. `IMBasisCollection` remains an exact page-to-basis mapping; it should not become a heterogeneous wave/APV/boundary/neutral state container.

The capable family protocol does not imply identical mathematics for every family. `IMHydrostaticModesBasis.discreteTransform(grid)` supports fixed-grid rigid-lid hydrostatic F/G. Native basis projection recipes retain their scalar pairings, and neutral labels retain the prescribed APV projector. Full-state waves and endpoint families can require several variables or derivative traces. Being able to evaluate F and G does not establish that either scalar component alone determines every physical coefficient.

## Four forward operations with different meanings

Write a sampled modal state as

$$y = A a,$$

where the columns of synthesis matrix $$A$$ are observations of the continuous modes, and $$a$$ retains the same scientific labels and normalization on every grid. Let $$B$$ discretize the continuous physical coefficient functional. It includes the required normalization, physical factors, volume pairings and endpoint terms. Let $$R$$ denote the selected recovery matrix. These are different objects.

| Proposed policy name | Definition | Meaning |
| --- | --- | --- |
| `physicalDual` | $$a_{\mathrm{physical}} = B y$$ | Approximate the manuscript's physical projection. Preserve quadrature/source errors; do not claim an exact retained-sample inverse. |
| `sampledLeastSquares` | Recover from the positive-observation-metric weighted synthesis using QR/SVD | Fit observed state samples. It is an exact left inverse on an observable retained span, within numerical accuracy. Its off-span action is determined by the declared observation norm. |
| `correctedPhysicalDual` | $$R = (B A)^{-1} B$$, constructed with stable left preconditioning | Recover through the physical test space. This is a Petrov–Galerkin recovery operator. Preserve $$B$$ separately; the correction changes the action on general inputs. |
| `truncatedCollocation` / extended-span recovery | Recover a declared larger span, then retain selected coefficient rows | Reject omitted columns that are explicitly represented in that larger span. The extension, retained labels and invertibility must be explicit. |

These names describe a recommended contract, not newly implemented APIs. Existing `forwardOperator="galerkin"` in the positive hydrostatic channels is weighted least squares. Existing root `"truncatedCollocation"` is the extended-span construction. A generic coupled overlap $$B A$$ can be complex and nonsymmetric; it is not a physical Gram matrix.

Introduce fixed coefficient amplitude scales $$C$$ and an observation whitening map $$D$$ so that $$a = C \widetilde a$$ and $$\widetilde y = D y$$. Then

$$\widetilde A = D A C, \qquad \widetilde B = C^{-1} B D^{-1}.$$

Rank and noise tests must declare these scales. For full-column-rank $$\widetilde A$$, let $$Q = \widetilde A \widetilde A^{+}$$. Corrected weak recovery satisfies

$$ (\widetilde B\widetilde A)^{-1}\widetilde B = \widetilde A^{+} + (\widetilde B\widetilde A)^{-1}\widetilde B(I - Q). $$

Thus exact recovery operators agree on the retained sampled span and differ on unresolved input, noise and nonlinear products. Among exact left inverses, least squares minimizes sample-noise amplification for the specified positive observation norm; multiplying its covariance comparison by any fixed output seminorm factor preserves that ordering. This does not make least squares the continuous physical projector or the best operator for every omitted-mode inventory.

State and source observations also have different contracts. A state layout may diagnose vertical velocity from horizontal divergence and continuity. An independent vertical-momentum source is not constrained by that state relation and needs its own source observation. WVM should assemble these physical polarizations and source recipes; the provider should accept generic prepared observations without depending on WVM classes. This follows the ownership in [AssessmentAPI.md](AssessmentAPI.md). Existing geostrophic state/source paths already distinguish corrected state recovery from raw physical source pairings.

## A correction to the earlier overlap diagnosis

The [DD grid requalification](NeutralLabelGridRequalification.md) found that the exact 65-point default rejected the common `[16 8 16]` physical-family allocation, despite full-rank state samples. The present investigation identifies the cause more precisely: numerical scaling of the weak-test equations produces a misleading rank rejection. This case does not demonstrate loss of retained information by the weak dual.

For its 69 physical plus two neutral amplitudes, the dimensionless sampled synthesis has rank 71 and condition number approximately $$3.80\times10^{6}$$. Coefficient amplitude normalization alone leaves the weak-test rows poorly scaled:

| Recovery equations for the same weak operator | Rank at unchanged relative cutoff `1e-12` | Condition number | Retained round-trip operator error |
| --- | ---: | ---: | ---: |
| Original overlap | 70 / 71 | `2.20e13` | `9.82e-6` if directly solved despite the rejection |
| Row-equilibrated overlap | 71 / 71 | `2.32e8` | `4.45e-8` |
| Orthonormal weak-test row frame, independently constructed by QR | 71 / 71 | `1.59e7` | `1.10e-9` |

For $$M = \widetilde B\widetilde A$$, invertible left preconditioning preserves the mathematical operator: $$(T M)^{-1}(T\widetilde B) = M^{-1}\widetilde B$$. No directions or physical observations were removed, and no acceptance tolerance was relaxed. Independently constructed QR and SVD row frames agree in the recovered map to relative `2.14e-15`. The row-scaled defining-equation residual is `2.11e-14`.

A second audit permuted the rows and multiplied them by complex phases and diagonal factors spanning twelve orders of magnitude. The original diagnostic rank changed from 70 to 41. Equilibration and row-frame recovery retained rank 71, preserved the recovered map to approximately `7e-13` and `1.4e-14`, and preserved the smallest trial/test principal cosine, `0.00508804665`. Future construction should equilibrate rows before the row-space rank decision as well; thresholding a raw poorly scaled dual would recreate the problem.

The preconditioned weak operator passes the independent refined-reference, physical error, noise, mandatory boundary and kernel-budget gates for this fixed allocation. Its fresh supported energy/enstrophy errors are approximately `3.63e-9` / `2.76e-12`. The historical maximum-count search remains a result of its original operator and rank policy. This investigation did not repeat the global count search with the improved construction.

## Recovery policy has measurable consequences

The coupled experiment held each grid, count allocation, continuous basis, physical/source scales and acceptance budgets fixed. It compared direct physical duals, the existing weak correction, preconditioned versions of that same weak correction, and sampled QR/SVD inverses. Construction and independently refined physical references were generated afresh. The source problem is the 1000 m neutral-cap example documented in the requalification report.

On the exact default grid, the following operators pass retained-recovery gates. The noise columns use the study's declared positive observation norm and separate physical energy/enstrophy output seminorms. Omitted-input columns use all available independently refined, nonmandatory labels outside each retained set. They are diagnostic responses, not accepted projection-error bounds.

| Physical allocation `[internal mixed APV]` | Recovery | Energy / enstrophy noise | Omitted-input energy / enstrophy response |
| --- | --- | ---: | ---: |
| `[16 8 16]` | Preconditioned physical dual | 179.11 / 224.77 | 115.23 / 0.2804 |
| `[16 8 16]` | Sampled QR | 179.10 / 11.35 | 114.33 / 0.9354 |
| `[41 5 37]` | Preconditioned physical dual | 9.9045 / 112.30 | 2.1841 / 0.7419 |
| `[41 5 37]` | Sampled QR | 9.9045 / 21.642 | 2.9769 / 0.8357 |
| `[46 6 42]` | Preconditioned physical dual | 266.71 / 632.66 | 268.50 / 21.311 |
| `[46 6 42]` | Sampled QR | 266.71 / 563.03 | 268.50 / 27.001 |

Least squares lowers the enstrophy noise gain, but the weak operator can better reject this omitted-input inventory. There is no universal winner. The common allocation is not a subset of the maximum allocation: it contains eight mixed-layer structures rather than six. Compare methods within a table allocation; comparisons between allocations also change the omitted-source set. The separate requalification report provides the identical-tail comparison across reduced bands.

The compact experiment covered five fixed grids and three allocations, not a new exhaustive search. Preconditioning does not make an unsuitable grid pass. The raw bandwidth grid's common allocation still fails physical boundary resolution and recovery/noise gates; the sampled QR energy noise gain is approximately `1.62e10`. This is distinct from the recoverable equation-scaling defect on the default grid.

## Root and non-root grids fit the same contract

A separate control used one converged common hydrostatic basis on five grids: root, physical depth, smoothed bandwidth, WKB-native and an explicit custom copy of the root points/weights. It retained 12 internal modes plus barotropic F at 33 points, with exponential positive stratification over 4 km. Independent analytical exponential modes supplied retained and omitted inputs. The common solve was independently refined from 256 to 384 spectral points: maximum shape error `9.30e-12`, equivalent-depth error `8.82e-12` and normalized root displacement `4.89e-14`.

| Grid / forward operator | Maximum F/G Gram distortion | Analytical retained recovery error | Response to omitted `j=13:31`, F / G |
| --- | ---: | ---: | ---: |
| Root / Galerkin | `4.87e-6` | `2.24e-9` | 0.975 / 0.991 |
| Physical depth / Galerkin | 0.0369 | `2.24e-9` | 1.358 / 1.389 |
| Smoothed bandwidth / Galerkin | `4.23e-5` | `2.24e-9` | 0.636 / 0.639 |
| WKB-native / Galerkin | `4.39e-5` | `2.24e-9` | 0.628 / 0.625 |
| Root / truncated collocation | `4.87e-6` | `2.24e-9` | `1.46e-12` / `2.14e-13` |

All Galerkin construction round trips were below `1.1e-15`; truncated collocation was below `5.3e-15`. Compatible F/G fields recovered the same internal coefficients within `8.3e-16`. The inactive barotropic G row was exactly zero. The custom copy produced bitwise-identical Galerkin operators: the grid's `kind` label does not change fixed-matrix numerics.

The two root policies use exactly the same continuous basis, retained sampled columns, points and weights. Truncated collocation suppresses the explicitly represented omitted span; this is a useful alias-control policy, not a different grid. On a higher folding band, `j=52:64`, its response remains approximately 1.413 F / 1.007 G, essentially equal to root Galerkin. It is not universally alias-free.

Physical-depth sampling illustrates a second distinction: it recovers the retained analytical coefficients accurately while its quadrature has 3.69% Gram distortion. Round-trip success cannot certify physical energy quadrature. These are bounded controls for one positive-stratification problem, not a universal grid ranking or timing benchmark. Their separately archived sample-noise gains use positive control-volume scalar F/G norms, not the coupled experiment's full-state energy/enstrophy norms.

## Positive recovery metrics and physical pairings must remain separate

There is also a smaller production numerical issue to fix. `IMProjection` currently diagnoses sampled-Gram rank before entering its positive diagonal-metric QR path. For

$$\Phi = \begin{pmatrix}1 & 1 \\ 0 & 10^{-8}\end{pmatrix}, \qquad W = I,$$

the sampled matrix has rank two and condition number `2e8`, but forming $$\Phi^{*}\Phi$$ rounds away the `1e-16` increment and produces rank one. The current constructor therefore takes the deficient-Gram branch and has round-trip error one; direct sample QR recovery has zero round-trip error in this example. A QR path entered only after a full-rank Gram test cannot prevent this failure.

For positive observation metrics, rank and recovery should be decided from the whitened synthesis directly, before forming a Gram matrix. Retain Gram distortion as a separate physical accuracy measurement. The example diagnoses a recovery path; it is not a claim that this example passes continuous-quadrature qualification.

A signed physical metric is different: its pairing can genuinely be degenerate even with full-rank samples. Do not silently substitute a positive metric, take absolute physical energies as norms, or reinterpret a signed Gram failure as a sample-rank failure. Neutral-label directions have zero physical energy, so physical energy alone cannot define their observation norm. The declared observations and positive norm must resolve these directions, for example through displacement or, at nonzero rotation, APV.

Likewise, volume quadrature weights are not automatically the right noise weights. Root rules can assign zero volume weight to endpoints whose traces remain essential observations. Record positive observation weights/scales independently; do not silently delete endpoint observations or add an arbitrary floor to a physical quadrature rule.

## Grid calculus and WVM integration

The same fixed-grid recovery engine can serve every grid, but a coherent model also needs a declared calculus. Points and quadrature alone do not define interpolation, differentiation, antiderivatives or one-sided traces at neutral joins. A local polynomial derivative, a native spectral quadrature and a separately reconstructed weak pairing do not automatically obey the same discrete integration-by-parts identity. Record the representation and test the required adjoint/boundary identities. Do not promise an exact SBP identity unless the actual operators satisfy it.

For WVM v4, the root path is structurally appropriate: inactive displacement boundaries leave interior G observations, while F includes the barotropic column. `verticalProjectionOperatorsWithRigidLid` already inverts a larger root-grid span then truncates and keeps physical moment weights separately. The new common basis can supply those modes without preserving every old numerical detail. Integration still needs explicit barotropic activity masks and analytic special-sector calculus: v4 uses a finite placeholder `h(1)=1`, whereas the common basis correctly records `h0=Inf`. Multiplying by h before excluding zero G can otherwise introduce `Inf*0`.

For active-boundary v5, several families must share one physical grid and keep independent boundary responses. The smoothed-bandwidth grid remains a reasonable provisional design after DD requalification, but it is not a complete transform or calculus contract. Current WVM free-surface construction still assumes WKB-native mapped-Chebyshev nodes for differentiation, interpolation, metric operators and its polynomial nonlinear-bandwidth heuristic; the balanced-state constructor explicitly rejects a custom scientific grid. Substituting new points alone would break those assumptions. Current Boussinesq construction also retains its positive-stratification restrictions; provider neutral-mode support does not remove them.

## Recommended implementation sequence

1. **Create one canonical recovery construction path.** Evolve existing `IMDiscreteTransform` into the prepared sampled-transform container, with an explicit named path such as `fromObservationMatrices` for arbitrary finite complex observation matrices, unique string column labels, sample/column counts, activity masks, observation metadata and declared scaling. Keep its existing scalar constructor behavior and `inverseMatrix` synthesis alias. Use a shared stateless numerical kernel for QR/SVD and preconditioned overlap construction. Keep physical dual, physical coefficient forms, overlap and recovery matrices separately named. `IMProjection` retains physical pairing/Gram semantics; `IMInternalModesDiscreteTransform` remains the aligned F/G adapter. Generic observation rows need not correspond one-to-one with z points. Existing scalar Gram/quadrature properties require explicit capability semantics on the generic path; an unavailable physical form must report unsupported rather than an apparent zero error.
2. **Fix the demonstrated numerical defects within that path.** Make positive-metric sample rank independent of Gram formation. Equilibrate physical test rows before rank decisions, orthonormalize their row space as needed and solve the resulting overlap without changing the physical dual. Measure sample conditioning, row-space/trial-space stability, round trips and defining-equation residuals separately. Test invariance under units, column normalization, row rescaling and ordering. A truly deficient declared test space must still fail explicitly.
3. **Use identical policies and assessments across grid designs.** Grid factories supply points/weights/provenance; fixed-grid construction performs no new EVP solve or hidden weight fit. Record any extended recovery span separately from the retained span. Keep source projection and state recovery as explicit operations. Do not optimize a grid using retained round trips alone or silently raise accepted counts after changing the numerical rank policy.
4. **Qualify the family and observation matrix.** Cover scalar hydrostatic F/G, full-state signed wave polarizations, inertial sectors, APV, zero-APV boundary responses and DD neutral labels. Exercise positive, neutral and signed cases; active/inactive boundaries; both rotation/frequency signs; empty sectors; repeated pages; real/complex data; independent source channels. Require independent retained-mode checks, exact inactive masks, boundary-shape checks and separate energy/enstrophy kernel budgets. Assess identical omitted inputs and representative nonlinear products in addition to retained recovery.
5. **Migrate WVM with compatible calculus and measured costs.** First compare root construction against the v4 contract, including barotropic integrals and physical source moments. Then replace v5's implicit native-grid assumptions with explicit reconstruction/derivative/trace capabilities and qualify active-boundary identities and products. Release/export the provider before consumer dependency adoption. Measure setup, application time and memory separately; this investigation does not establish a fastest grid or highest usable nonlinear resolution.

The immediate recommendation is the common recovery engine and the two demonstrated rank/conditioning fixes. Use sampled QR/SVD as the explicit observability/noise baseline, keep the manuscript source dual intact, and retain preconditioned physical-dual and extended-span recovery as distinct policies. Only after that comparison should a production recovery policy and retained bandwidth be selected for each WVM use case.

The proposed `IMDiscreteTransform` expansion needs an explicit capability contract. Generic forward/back dimensions use `sampleCount` and `columnCount`; string `columnLabels` are canonical, while z, quadrature weights and numeric `modeNumber` become optional scalar-adapter metadata. Physical pairing diagnostics belong in an optional record: absent metrics/Grams remain unavailable, Gram errors/ranks/conditioning report unsupported, and no placeholder positive-definiteness claim is made. Keep sample rank, test-space rank and overlap conditioning separately named. If preserving the legacy scalar property contract prevents those semantics, use a small generic prepared container behind the existing adapters instead. Reusing one construction engine matters more than forcing unlike mathematical objects into one class.

## Evidence, reproduction and verification ledger

The relevant production paths are the [fixed-grid hydrostatic basis bridge](../@IMHydrostaticModesBasis/IMHydrostaticModesBasis.m), [qualified root factory](../@IMInternalModesDiscreteTransform/hydrostaticRoot.m), [grid constructors](../@IMVerticalGrid/IMVerticalGrid.m), [positive-metric projection constructor](../@IMProjection/IMProjection.m) and [current Gram-based rank initialization](../@IMProjection/initializeProjection.m). Consumer checks used [v4 rigid-lid operators](../../wave-vortex-model/@WVStratification/verticalProjectionOperatorsWithRigidLid.m), [v5 balanced-state setup](../../wave-vortex-model/+WVInternal/buildFreeSurfaceBalancedState.m), [v5 wave-state setup](../../wave-vortex-model/+WVInternal/buildFreeSurfaceWaveState.m) and [native-grid calculus](../../wave-vortex-model/+WVInternal/qgVerticalOperators.m).

Local research archive: `literature/wave-vortex-model-coordinate-grid/data/transform-investigation-20261102/`. It remains ignored by Git, like the preceding grid archives. It contains the exact executed MATLAB probes, CSV tables, saved matrices, the source study configuration, provider/study/WVM revisions and working patches, current MATLAB source hashes, manuscript hash and a SHA-256 file manifest.

| File | Evidence |
| --- | --- |
| `coupled-operators.csv`, `coupled-conditioning.csv`, `coupled-results.mat` | Five fixed grids, three physical allocations, six operator constructions, independently refined validation, sample/noise/tail diagnostics and saved decisive matrices. |
| `investigateTransformOperators.m`, `study-configuration.mat` | Exact coupled experiment and construction/refined-reference settings. Uses the prior completed 65-point comparison for fixed grid geometry; mode references are freshly solved. |
| `audit_rows.m`, `audit-rows.csv`, `audit-rows.mat` | Independent QR row-frame derivation, left-preconditioning and equation-residual checks, row-scale/phase/permutation invariance. |
| `positive-scalar.json` | Positive-metric sample-rank versus rounded-Gram control. |
| `root_grid_coherence_probe.m`, `.csv`, `.mat` | Same-basis grid/policy comparison, analytical retained/tail inputs, custom-grid equality and independently refined root-solve evidence. |
| `provenance.json`, `sha256.json` | Source/revision context and exact artifact checksums. |

The scripts preserve their executed workspace and temporary-output paths for provenance. To rerun elsewhere, change those setup paths. The coupled probe needs the prior `Nz65/configuration.mat` and `comparison.mat` documented in [NeutralLabelGridRequalification.md](NeutralLabelGridRequalification.md); the configuration is also copied here. The row audit consumes `coupled-results.mat`. The analytical root probe is independent of the neutral-cap archive. Use the workspace's provider dependencies and `matlab -batch`; on this local host MATLAB must run outside the Codex sandbox. One computational thread was used for the coupled experiment; this is setup, not a different mathematical policy.

Verification: both numerical experiments and the independent audit completed successfully. Code Analyzer reported only finite-loop dynamic-array growth advisories in the two exploratory probes (12 and two messages); the independent audit had zero messages. No production MATLAB source, package manifests, dependency declarations, released snapshots or generated website output changed in this investigation. No production regression suite was rerun for the report-only change. Markdown link/structure, source scope and whitespace checks were performed. The repository has no configured `docs:check` task; this standalone research note does not feed generated website documentation. This is not a full family qualification, a new maximum-band search, a performance benchmark or completed WVM integration.
