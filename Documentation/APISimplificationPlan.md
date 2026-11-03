# V2 API integration and simplification plan

Status: implemented and verified on 2026-11-03. PR [#68](https://github.com/JeffreyEarly/internal-modes/pull/68) now integrates the newer target-branch work with the capable-basis and manuscript endpoint APIs below. See the [integration verification record](APIIntegrationVerification.md) for the full-run/focused-rerun ledger and unchanged figure results. The [earlier verification record](APISimplificationVerification.md) describes the preceding consolidation; its historical checks are separate.

## Goal and scope

Provide one public workflow: `IMBackground` → scientific family problem → solver → capable family basis. Integrate the useful recovery, identity, provenance and normalization work added to the target branch without restoring `IMModes`, `IMModeBasis` or `IMModeBackground`. Remove `IMVerticalSampling`. Replace the general public endpoint-rotation framework with the two boundary-buoyancy basis choices used in the manuscript.

The provider authoring repository is `internal-modes-evp`, published as `JeffreyEarly/internal-modes`. PR #68 uses `v2/simplify-mode-api` at `f144924d72b1250f053e1886201e7c4880616926`, targeting `InternalModesEVP` at `89ab5a7501189d1bf11adb471a03b03e9781719a`. Their merge base is `8717aca878fff72b5e2c19f962e7dc6c8c435888`. The assessed merge has 23 conflicting files: 12 runtime files, five tests and six documentation files. Recheck these references before implementation if either branch advances.

Keep the current manuscript and free-surface WaveVortexModel integrations working. A complete port of `WVTransformBoussinesq`, a new sampling abstraction and package release/export work are separate projects. Do not edit released package snapshots under `OceanKit/` or publish website output during this integration.

## Public API after integration

### Construction, counts and evaluation

Keep direct family construction and capable numerical, analytical and specialized results:

```matlab
background = IMBackground(N2=N2,zDomain=[-D 0],rotationRate=rotationRate,latitude=latitude,surfaceCondition="freeSurface");
solver = IMSolverSpectral(nEVP=512,coordinateKind="bandwidth");
problem = IMInternalModes.geostrophicAPVModes(background,g0=g0,gd=gd);
basis = solver.solveEVP(problem,nModes=4);
F = basis.F(z,normalization="depthRMSF");
Gz = basis.evaluate(z,variable="G",derivativeOrder=1,normalization="depthRMSF");
```

- `nModes` is required and means exactly the total retained columns, including special sectors. Zero returns an empty basis. An unavailable exact inventory fails explicitly; never clamp the request or append implicit columns.
- Explicit branch requests use explicit matching count vectors. Preserve defined ordering. Translate older internal/positive counts by scientific inventory, rather than mechanically renaming options.
- Positive bulk pages require `nModes`; inertial pages use `nInertialModes`. A zero-only request needs only its inertial count. Preserve repeated wavenumbers and page order.
- Scalar, collection and chunked evaluation share named `variable`, `derivativeOrder` and optional `normalization` arguments. Require the variable; default derivative order to zero. Use F, G, APV or generic u as supported. Keep convenience methods such as `Gz(z)` where supported; do not revive derivative names as evaluator variables.
- Omitted columns mean all columns; explicit `columns=[]` means no columns. Preserve order and repeats where supported, and keep any distinct-column restriction of `selectColumns` explicit. Resolve the all-column sentinel before dispatch to specialized evaluators.
- Validate every selected collection page before invoking streaming callbacks. Empty numerical, analytical and specialized results preserve sample-by-zero shapes and descriptive metadata.

### Boundary-buoyancy basis choices

Keep explicit canonical construction and expose one scientific basis-changing method:

```matlab
problem = IMGeostrophicZeroAPVModes.atWavenumber(background,k,endpoints=["surface","bottom"]);
boundary = solver.solveGeostrophicZeroAPVModes(problem);
orthogonal = boundary.orthogonalizeGeneralizedEnergy(g0=g0,gd=gd);
hbb = orthogonal.boundaryDepth;
M = orthogonal.boundaryResponseMatrix;
```

The canonical basis has unit buoyancy response at each requested endpoint. `endpoints` chooses the inventory; `g0` and `gd` specify the generalized-energy pairing independently. A single requested endpoint is a complete one-dimensional family. For a free surface, the surface response is `B(0)=G(0)-F(0)`; the bottom response is `B(-D)=G(-D)`.

Preserve `endpoints=strings(1,0)` as a complete empty family with empty profiles, maps, forms and boundary depths. Orthogonalization of that originally empty family returns an empty result without an eigensolve. Selecting no columns from a nonempty family is a proper reduction and cannot subsequently be orthogonalized.

Require finite signed coefficients for requested endpoints. Accept positive infinity only for an absent endpoint, explicitly omitting its term before arithmetic; reject NaN and negative infinity. An infinite coefficient must never remove an endpoint or produce `Inf*0`. Keep composed-transform handling of inactive endpoints consistent with this rule.

`orthogonalizeGeneralizedEnergy` constructs the manuscript's generalized-energy orthogonal combinations in the complete originally chosen endpoint space. Its default scaling is the manuscript's unit endpoint-response geometry, rather than unit generalized energy. Require explicit `g0` and `gd`, consistently with the corresponding PV pairing. The operation performs no new boundary-value solve.

The recommended order is: choose endpoints → solve canonical responses → orthogonalize if needed → select columns → normalize → evaluate or prepare recovery. Orthogonalization of a proper reduced selection fails clearly. It must neither restore discarded directions nor report a restricted Rayleigh quotient as a full manuscript boundary depth. A full-span prior normalization must not change the physical orthogonal directions.

Expose `boundaryResponseMatrix` as a read-only map with canonical-endpoint rows, current-mode columns and wavenumber pages:

$$a_{\mathrm{boundary}} = M a_{\mathrm{current}}.$$

For the canonical basis, M is the identity. For the default orthogonal basis, it is the manuscript coefficient matrix C. Normalization rescales its columns; selection can make it rectangular. Preserve the underlying canonical identity so compatibility does not depend on a display name.

Expose signed `boundaryDepth` for the orthogonal basis as the manuscript quantity, distinct from ordinary EVP equivalent depth `h` or eigenvalues. Preserve this intrinsic label under subsequent column scaling and selection; transform the actual current-coordinate Gram matrices separately. Canonical boundary responses have no ordinary EVP spectrum, and no boundary-depth label is implied before orthogonalization.

Keep physical-form getters in current coordinates: `energyMatrix`, `generalizedEnergyMatrix(g0=...,gd=...)`, `surfaceBuoyancyMatrix` and `bottomBuoyancyMatrix`. Do not introduce a public hierarchy of form objects. Apply the same coordinate map to F, G, B, derivatives and endpoint traces, and congruence transformations to all forms. Coefficients, physical duals and recovery operators must transform consistently.

Remove public `rotateWithPencil`, `rotateBoundaryDepth`, `rotateSurfaceBuoyancy`, `rotationRightMetric`, the general rotation registry, and the public `rotationName`/`rotationMatrix` contract replaced by the scientific operation and response map. Replace the older `h0` endpoint-depth name with `boundaryDepth`. Keep only internal coordinate algebra needed by the retained operations. Use manuscript terminology in public help and examples.

### Endpoint numerics

For nonzero horizontal wavenumber, form the real symmetric canonical generalized-energy matrix and solve

$$H_g C = C\,\operatorname{diag}\!\left(\frac{h_{\mathrm{bb}}}{2k^2}\right), \qquad C^T C = I.$$

Use a symmetric eigensolver in canonical endpoint coordinates, or an equivalent formulation with positive endpoint geometry. Generalized energy can be indefinite or singular: negative, zero and repeated boundary depths are valid. Do not require a Cholesky factorization of H_g or manufacture positive signatures. Preserve the family's existing wavenumber admissibility rules.

Qualify symmetry, eigen-residuals, coordinate rank and degeneracy with scale-relative criteria. Compare invariant subspaces for repeated eigenvalues; do not require individual eigenvectors to remain identical. Fix ordering and signs deterministically where the spectrum permits, while recognizing that repeated eigenspaces have no unique directions. Physical energy need not diagonalize with generalized energy.

A zero generalized-energy direction remains a usable basis vector. A requested unit-generalized-energy normalization or inverse-norm projection must reject that null direction explicitly, without deleting it or replacing a zero divisor by one.

The scientific source is the manuscript's boundary-buoyancy construction, signed boundary-depth normalization and coefficient-conversion equations in `literature/ape-apv-free-surface/main.tex`. No manuscript equation change is required by this API cleanup.

## Work to retain from the target branch

### Capabilities, scientific identity and provenance

Provide one authoritative structured `capabilities` contract on every capable family, including generic scalar u bases. Describe supported variables, derivative orders, spectral quantities, normalizations, endpoint trace conventions and conditional projection/recovery requirements. Update scalar, collection, chunked and product consumers to use it consistently. Capability inspection describes supported operations; it does not predict numerical rank, acceptable noise or recovery success.

Preserve separate `scientificSpecification` and `constructionRecord` information. The former identifies the requested and resolved scientific family, background, EVP/energy role, support, endpoints and horizontal meaning. The latter records the original solver request, configured backend, original counts/branches, diagnostics, eigenpair/source indices and subsequent selections or basis choices. Selection and normalization must not rewrite the original request or discard qualification evidence. Default numerical evidence is unverified; resolution, requested count or an analytical class name does not certify convergence.

Distinguish fixed horizontal wavenumber, signed per-column `kSquared` at fixed frequency, horizontal independence and unspecified developer-supplied identity. Preserve evanescent roots and zero values where the family supports them. Real inertial bases have `k=0`; horizontally independent families must not acquire a fictitious zero wavenumber. Collections reject one scalar kappa for per-column dispersion, permit an explicit kappa for independent structures, and preserve repeated-page ordering.

Inspection of records, capabilities or available spectra must remain cheap and must not trigger normalization integrals or another solve.

### Normalization and selection

Keep common `withNormalization`, per-call overrides and custom rules on numerical, analytical and specialized bases. Preserve independent defaults on returned copies; route retained active-normalization setters through the same validation and cache hooks. Normalization is lazy and reuses existing modes and quadrature infrastructure.

Reject zero or nonfinite norm magnitudes and nonpositive or nonfinite custom scale factors. Signed scientific self-pairings use a positive square-root-of-absolute-value divisor and a separate sign. Define `normSignatures` as signs of the active normalization form; positive profile norms have positive signatures. Preserve signed physical forms independently of the chosen profile normalization.

Endpoint bases reserve `"endpointResponse"` and `"generalizedEnergy"` so custom registration cannot silently change the promised canonical geometry or signed-energy meaning. Other profile rules and custom names can still be replaced.

Selection must retain the original integration effort needed for high modes, while assessing normalization admissibility only on the retained columns. A discarded null barotropic G column must not prevent valid normalization of selected internal columns. Keep original-inventory quadrature metadata independently of the sampling class being removed.

Selection, custom-rule replacement and coordinate changes must invalidate or detach affected derived caches. Reuse unchanged rules and avoid additional solves, grid construction or weight fitting. Preserve prescribed neutral APV duals with the correct reciprocal normalization scaling, including explicit empty selections. Prepared transforms retain their captured basis; creating another normalization must not silently mutate an existing transform.

### Physical projections and numerical recovery

Retain `prepareIMRecovery`, `IMDiscreteTransform`, `IMScalarDiscreteTransform`, `IMProjection` and the useful aligned F/G recovery changes. Keep physical pairings, continuous signed forms and positive majorants, physical source duals, and numerically corrected recovery operators distinct. Preserve `projectPhysical`, `prepareRecovery=false` and recovery-failure diagnostics when a physical dual remains usable.

Use strict rank qualification and an explicit fitting/noise model. Determine rank from the observation factorization, such as QR or SVD, rather than rounded Gram eigenvalues. Do not label a pseudoinverse fallback as successful full recovery. Preserve caller-supplied points and weights exactly.

Move `projectionOnGrid` and `scalarTransform` from the removed wrapper onto capable families through shared internal helpers. Require an explicit variable and call the named evaluator. These adapters perform no eigensolve or weight fitting. Retain physical scalar recipes only where the source and metric information supports them: hydrostatic F/G, geostrophic APV F/G when their required observations are supplied, MDA G and neutral APV. Keep unsupported MDA F, wave/inertial/fixed-frequency, zero-APV, generalized-potential-enstrophy, derivative and companion-trace physical recipes explicit. An auxiliary EVP scalar pairing for wave G is not a full ocean-state physical projector; generic caller-assembled matrix recovery remains available independently.

For a selected physical subspace, assemble and solve its retained target-form system. Do not obtain its physical dual by merely slicing rows from a full-space dual; those maps generally differ for nonorthogonal forms. Preserve independent endpoint observations and their literal coefficients, including zero, infinite and constrained limiting cases. Never infer an arbitrary input's missing endpoint or companion trace from a retained mode's boundary law.

### Root transforms and composed geostrophic transforms

Repair the semantic mismatch in `hydrostaticRoot`: the returned capable hydrostatic basis must accept the prepared recovery and forward-operator injection used by the root constructor. Port that bridge into native hydrostatic and analytical fixed-grid paths without restoring the wrapper. Preserve full extended-span preparation followed by output selection, prefixing, reweighting and qualification provenance. Keep retained-span Galerkin and full-span truncated-collocation policies distinct and explicit.

Use scientific column activity, combined with actual labels, in projection recipes. Hydrostatic barotropic G is inactive even if source-EVP form metadata exists; MDA's infinite-depth constant G is active. Do not infer activity from sample magnitudes or produce `Inf*0` for an inactive column.

Migrate `IMGeostrophicTransform` from `rotationName`/`rotationMatrix` checks and conversions to canonical scientific identity and `boundaryResponseMatrix`. Keep its full-active-endpoint requirement, consistent PV/endpoint generalized-energy coefficients, and source-versus-state distinctions. For an ocean-state projection, subtract the reconstructed PV endpoint anomaly before determining boundary-buoyancy coefficients. Independent forcing/source observations remain independent inputs. Verify that changing endpoint coordinates preserves reconstructed states, physical energy cross terms and source/state projectors, including complex fields and trailing dimensions.

## Complete removal of prescribed sampling

Remove the incoming `@IMVerticalSampling` class and evaluator, `ExamplesV2/PrescribedVerticalSampling.m`, `UnitTestsV2/IMVerticalSamplingTests.m`, `UnitTestsV2/Fixtures/IMPrescribedStaggeredSampling.csv`, and its fixture README entry. Remove `Documentation/PrescribedVerticalSampling.md` and `Documentation/PrescribedVerticalSamplingVerification.md` and their links or claims in documentation/example indexes, `CommonScalarProjection.md`, `SharedRecovery.md`, current changelog entries and class comments.

Search for registrations and generated references. Remove sampler-only `samplingIntegrationInfo` and `samplingBreakpoints` declarations/helpers when no retained consumer requires them; preserve normalization integration budgets through an independent internal mechanism.

Point sampling is ordinary `basis.F(z)`, `basis.G(z)` or `basis.evaluate(...)`. Retain `IMVerticalGrid`, caller-assembled observation matrices, the common recovery engine and generic observation-recovery tests. Do not add a replacement observation-layout class, cell-average abstraction or standalone sampling helper in this scope.

## Scientific invariants and consumer scope

- All public wave source paths retain `F = D Gz`, including fixed wavenumber, fixed frequency, inertial and exact zero-depth paths. Do not copy the incoming wrapper conversion onto already-converted native profiles. Keep the hydrostatic and balanced defining relations distinct.
- WaveVortexModel retains its single provider-to-consumer conversion `F_WVM = (h/D) F_provider`, explicit wave normalization and `"kConstant"` use. Preserve reconstruction, source pairing, calculus, quadrature, energy and restart behavior.
- Keep `rotationRate`/`latitude` support. `"wkb"` continues to mean the classic N coordinate; `"bandwidth"` remains the separate rotation-aware coordinate. No new nontraditional equations are introduced.
- Preserve positive Dirichlet-Dirichlet neutral sine labels, no constant neutral label, exact counts and explicit interval traces. Keep the neutral continuous APV source dual distinct from sampled coefficient recovery.
- Keep canonical endpoint construction in manuscript Figure 6 and the existing free-surface WVM builder. Neither consumer needs to adopt orthogonal coordinates merely because the operation exists.

The WVM source migration is already committed as `c4621b3e258be159fc9cff101203118a53a7f764` on local branch `feature/im-v2-api-migration`. Audit changes required by the integrated provider and keep further edits bounded to these existing consumers; do not begin the broader `WVTransformBoussinesq` port.

The manuscript source/figure migration is already committed as `82fb249` on `master`. Overleaf work stays on that default branch. Preserve its seven PDF baselines and plotted arrays for the final comparison. No Overleaf push is implied by this plan.

## Implementation sequence and acceptance criteria

### 1. Preserve state and establish the integration baseline

Record current branch tips, the merge base, conflict inventory and consumer commits. Preserve the user's current uncommitted free-surface example edit and this planning document separately. Use an isolated integration checkout if needed so the merge cannot overwrite that example. Do not reset or silently restore it.

Capture the target API examples above as acceptance cases. Reuse available numerical, figure, root-grid and consumer baselines; record inherited failures before calling them integration regressions. Merge the target into the existing PR branch, preserving both histories. Resolve by final contract, without a blanket choice of either side. Check automatically merged code for semantic breaks as well as marked conflicts.

Acceptance: every change is attributable to the earlier consolidation, retained incoming functionality or an explicit decision in this plan; the user's unrelated local edit is preserved.

### 2. Integrate the capable-basis foundation

Reconcile structured capabilities, scientific and horizontal identity, construction records, exact counts, named evaluation, explicit-empty behavior, normalization and selection. Make numerical, analytical and specialized bases satisfy the applicable contract. Adapt collection/product consumers and remove incoming dependencies on the retired namespace/wrapper.

Acceptance: consistent scalar/collection/chunk behavior, correct empty shapes and page order, no implicit extra columns, truthful capability rejection before callbacks, cheap inspection, lazy norms and no extra solves. High selected modes retain integration accuracy; discarded null columns cannot poison their normalization.

### 3. Implement the manuscript endpoint API

Implement canonical identity and the current response map, `orthogonalizeGeneralizedEnergy`, signed `boundaryDepth`, selection/normalization covariance and current-coordinate physical forms. Apply the same contract to supported numerical and analytical endpoint paths. Remove obsolete public rotation machinery and migrate composed geostrophic transforms.

Acceptance: canonical unit endpoint responses; coupled-equation residuals; orthonormal coefficient columns and diagonal generalized energy in default orthogonal coordinates; negative, zero and repeated depths handled correctly; complete empty and single-endpoint support; clear rejection after a proper reduction; explicit active/absent endpoint coefficient validation; no restored directions or new boundary solve. Coordinate changes preserve physical states and projectors.

### 4. Integrate recovery and transforms; delete prescribed sampling

Integrate shared recovery and native scalar adapters, repair prepared hydrostatic root/fixed-grid bridges, and reconcile scientific activity and projection recipes. Remove the complete prescribed-sampling feature and orphan support identified above. Retain generic matrix recovery.

Acceptance: unchanged physical source duals, qualified recovery failures, original point/weight use, correct root policies and extended-span behavior, inactive hydrostatic G0 versus active constant MDA G, and no production reference to any removed class or rotation method.

### 5. Migrate examples, tests and canonical documentation

Translate incoming tests by scientific intent. Keep useful capability, provenance, horizontal identity, normalization, scalar projection, observation recovery and endpoint coverage on capable bases. Use a test-only fixture where an old wrapper fixture supplied artificial developer data. Replace arbitrary public-rotation tests with manuscript invariants rather than preserving an unnecessary API to satisfy them. Delete only sampler-specific coverage alongside that removed feature.

Update executable examples, README, changelog, canonical inline API help and creation, basis-contract, normalization, endpoint-choice, projection and recovery guides with final before/after usage. Reconcile incoming documents such as `BasisContracts`, `BasisNormalization`, `BasisRotations`, `CommonScalarProjection` and `SharedRecovery`; describe supported scientific endpoint choices instead of a generic rotation framework. Preserve historical verification reports as historical evidence, and add an integration ledger for the new implementation.

Audit manuscript and WVM callers against the final API. Keep manuscript source edits on `master`, WVM edits on its migration branch and provider integration on PR #68's branch. Keep user-local example changes outside the integration unless separately reconciled.

Acceptance: examples use the final API; tests preserve mathematical coverage; active documentation contains no retired usage; historical reports do not claim current compatibility or newly completed verification.

### 6. Verify the coherent implementation

Start with affected focused suites, repairing and rerunning only failed methods. Include capability/provenance/horizontal identity, native protocol/counts, evaluation/collections, normalization/cache covariance, endpoint/geostrophic composition, recovery/observation/scalar projection, and hydrostatic root integration tests. Cover numerical, finite-difference, analytical and special families where applicable. Then run the complete provider `UnitTestsV2` suite after integrated public contracts are coherent.

Exercise explicitly: negative/zero/repeated generalized-energy directions; endpoint equations and maps; coefficient and form covariance; complex fields/trailing dimensions; selected high harmonics; discarded null columns; custom-rule replacement/copy independence; physical off-span source projection versus corrected recovery; explicit-empty requests; unsupported companion traces; both root-transform policies; and no-repeat-solve behavior.

Run affected examples once after the coherent batch. Run bounded existing WVM checks for construction, reconstruction, energy/source projection, forced evolution, restart, stored-mode reuse and convergence/shooting qualification. Repeat the archived physical-quantity comparison where changed production paths warrant it, and run the WVM artifact policy check. Full WVM qualification remains separate from these bounded checks.

Regenerate the seven manuscript figures once after the final coherent mathematical batch. Compare physical arrays and metadata first, then rendered figures against preserved baselines, and inspect the results. Report every scientific or visible change. An arbitrary direction within a repeated eigenspace is assessed by its subspace; it must not be mistaken for a changed physical projector. The goal for current manuscript figures is unchanged scientific content, but do not assume that outcome before comparing.

Run Code Analyzer once on final changed MATLAB files, then only on files affected by repairs. Generate affected canonical API documentation once into temporary output after the help-text batch; check Markdown fences, local links, removed symbols/callers, whitespace, repository scope, package manifests and generated-artifact policy. There is no configured provider `docs:check` task: record bounded checks actually run. Maintain a ledger and avoid repeating successful gates without a relevant subsequent change.

Acceptance: the full provider suite and required bounded consumer checks pass, or actual inherited failures/missing assets are named precisely. Figure differences, residual risks and unrun qualification are reported explicitly. Historical passing counts are not presented as a new integrated full-suite result.

### 7. Deliver the reviewed integration

Keep commits reviewable around coherent contracts and their callers; preserve merge history. When executing the authorized integration, commit provider changes, push the existing PR branch and update PR #68 to describe the final API and actual verification. Recheck the base and conflict status. Keep the PR draft until ready for review; no force-push, default-branch merge, new issue, WVM PR or release is required by this plan.

Commit any newly required bounded WVM or manuscript edits in their existing repositories and appropriate branches. Preserve the user's local edits and distinguish provider integration from downstream changes. Follow the shared external-communication policy for any tracker action. The user has already authorized a PR in `JeffreyEarly/internal-modes`; no additional tracker workflow is needed merely to formalize this plan.

Acceptance: one coherent provider PR, separately identifiable consumer changes, a verification ledger, an explicit figure-change report and no unexplained generated artifacts or package snapshot edits.

## Separate release follow-up

The current WVM dependency pin still refers to the older beta package lacking `IMBackground` and the final wave convention. Source-checkout integration alone does not make that exported package graph usable. After the provider API is reviewed and released, update the consumer dependency floor and CI pins through the public MATLAB MPM workflow, then qualify a clean exported-package installation. Do not hand-edit manifests, release/tag packages or update OceanKit snapshots during conflict integration.

## Planning verification ledger

- Reviewed shared package-design, documentation, caching and package-release guidance, together with applicable MATLAB/inline-documentation guidance from the earlier implementation.
- Inspected branch histories, conflict and automatic-merge results, native and incoming basis/recovery contracts, manuscript equations, endpoint consumers and bounded WVM integration paths.
- Reconciled independent read-only reviews of normalization/identity, recovery/consumers and caller/documentation/removal scope.
- Checked the final plan's Markdown fences, local links, conflict markers and whitespace. Preserved the preexisting free-surface example edit.
- This planning update changes only this Markdown document. No MATLAB source, tests, figures, package metadata or website output were changed, and no numerical verification was rerun.
- Historical baseline results remain in the earlier verification record. The pending integration needs the new verification sequence above.
