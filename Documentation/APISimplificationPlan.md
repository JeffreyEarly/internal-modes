# Simplify the V2 public API

Status: implemented locally on 2026-11-02. The five coordinated API changes, consumer migrations and local verification are complete. See the [verification record](APISimplificationVerification.md) for the full-run/focused-rerun ledger, unchanged figure comparison and bounded consumer checks. The staged scope and acceptance criteria below record the agreed plan.

## Agreed direction

Make the existing physical background → family problem → solver → capable family basis workflow canonical. Remove the public `IMModes` namespace entirely. Retire the parallel `IMModeBasis` wrapper after relocating its useful functionality into capable family bases and shared internal helpers. Keep specialized family classes where their mathematics differs; do not force endpoint responses, neutral labels and exact zero-depth kernels into a finite-eigenvalue scalar EVP.

The five coordinated changes are:

1. Make `nModes` count the total retained columns, including surface, negative, zero and barotropic sectors. Do not append implicit columns beyond a request. A branch-specific request must count its final selected inventory. Use a consistent explicit error when an exact requested inventory cannot be resolved.
2. Return capable family bases with their applicable diagnostics, derivatives, normalization, selection, projections, transforms and endpoint rotations. Convenience construction must not hide these operations behind a private source object.
3. Preserve per-evaluation normalization overrides, an active default normalization and custom rules. Compute expensive normalization factors only when needed and reuse them without another eigensolve. Selection and normalization must preserve aligned F/G fields, scientific identities, duals and metric signatures.
4. Select canonical zero-APV endpoints explicitly with `endpoints`. Generalized-energy coefficients specify a problem, pairing or basis rotation; they do not indirectly select canonical responses.
5. Define the manuscript wave convention, F = D Gz, in the wave family itself. Update the forward and inverse diagnostic relations, derivatives, normalization, endpoint terms, projection duals and transforms together. Every public construction path for that wave family must agree. Hydrostatic and geostrophic families retain their distinct defining relations.

Preserve the current `rotationRate`/`latitude` inputs, physical defaults, classical N-based `"wkb"` coordinate, separate `"bandwidth"` coordinate and support for full rotation information. This cleanup does not introduce additional nontraditional equations. Preserve the current fixed-endpoint neutral-label sine family and its exact zero-count behavior.

## Intended usage

The following illustrates the implemented canonical API.

```matlab
background = IMBackground(N2=N2,zDomain=[-D 0],rotationRate=rotationRate,latitude=latitude,surfaceCondition="freeSurface");
solver = IMSolverSpectral(nEVP=512,coordinateKind="bandwidth");
problem = IMInternalModes.geostrophicAPVModes(background,g0=g0,gd=gd);
basis = solver.solveEVP(problem,nModes=4);
F = basis.F(z,normalization="depthRMSF");
G = basis.G(z,normalization="depthRMSF");
gram = basis.gramMatrix(variable="F");
```

Canonical endpoint responses remain explicit:

```matlab
problem = IMGeostrophicZeroAPVModes.atWavenumber(background,k,endpoints="surface");
basis = solver.solveGeostrophicZeroAPVModes(problem);
```

Keep existing explicit physical-input constructors where useful. Route both forms through one physical problem implementation, with conflicting sources rejected. Analytical and numerical backends must return the same scientific family contract; an analytical marker that only works through the removed namespace is not sufficient.

The count and branch-selection signatures must be specified and tested before implementation spreads to callers. Avoid introducing multiple competing count policies or a new public selection-object hierarchy solely to replace the removed facade.

## Baseline and scope

The authoring checkout is `internal-modes-evp`, whose remote is `JeffreyEarly/internal-modes`. At planning time its branch is `v2/neutral-label-manuscript-boundaries` and HEAD is `8717aca878fff72b5e2c19f962e7dc6c8c435888`. It has existing uncommitted neutral-label, basis, collection, documentation and test changes. Preserve that working state as the scientific baseline and record inherited changes separately before any implementation edits. Do not reset the branch or accidentally drop those changes while creating an isolated checkout.

The manuscript repository `literature/ape-apv-free-surface` also has the preceding eight-script migration, seven regenerated PDFs and its comparison record as local changes. Preserve its numerical baseline and distinguish subsequent consolidation edits from that earlier migration.

Do not edit released package snapshots under `OceanKit/`. Do not publish website output, change package versions, push source, create tracker communications or publish a release as an incidental part of this local refactor. Package metadata changes, if required, must use the public MATLAB MPM API and the normal release workflow.

## Dependency inventory

| Area | Existing coupling | Required migration |
| --- | --- | --- |
| `@IMModes/IMModes.m` | Family construction, special-sector selection, analytical dispatch, exact zero-depth evaluator | Move scientific functionality to the family/problem/backend that owns it; then remove the namespace |
| `@IMModeBasis/` | Profile conversion, eager normalization, metadata, selection and hydrostatic sampling bridge | Move useful operations to capable family bases and shared internal helpers; then remove the wrapper |
| `IMModeSelection.m`, `@IMSolver/IMSolver.m` | Separate branch-count and automatic-special policies | Establish one total-count contract with explicit branch selection and exact failure behavior |
| `@IMSolverAnalytical/` | Marker interpreted only by `IMModes` | Supply a real backend path or use the existing capable analytical solution entry points consistently |
| `@IMBasisCollection/` | Type-specific facade handling and incomplete F-derivative metadata | Consume the capable family protocol and preserve supported evaluation/projection operations |
| `@IMInternalModesDiscreteTransform/hydrostaticRoot.m` | Calls the facade, returns a typed wrapper, uses its selection/root/transform methods | Preserve the qualified root transform and reusable returned continuous basis through the capable hydrostatic family |
| `@IMInternalModesDiscreteTransform/IMInternalModesDiscreteTransform.m` | Constructor enforces wrapper type | Validate alignment, normalization and supported basis capabilities without requiring the retired type |
| Exact zero-depth waves | Scientific construction exists inside `IMModes` | Add an explicit family construction path before deletion; preserve critical intervals, finite auxiliary count and manuscript convention |
| Neutral labels and endpoint responses | Existing capable specialized bases plus facade adapters | Preserve specialized construction, fixed-endpoint labels, prescribed duals, signed metrics and endpoint rotations |
| Manuscript scripts | Eight mode scripts, including a direct frequency-squared pencil and a nonlocal endpoint experiment | Use canonical construction without changing those experimental equations; rerun numerical and visual comparisons |
| WaveVortexModel | Uses lower-level wave profiles, bulk solves, source projections and restart data | Audit and adapt required dimensional conversions; verify physical reconstruction and projections rather than relying on normalized plots |

The provider also has facade references in `ExamplesV2/HydrostaticRootTransformComparison.m`, README, creation documentation and five focused test classes. Existing wrapper tests must be translated into scientific contract tests; do not delete their coverage along with the class.

## Delivery sequence

### Stage 0: freeze contracts and capture baselines

- Preserve tracked and relevant untracked working changes without altering the scientific baseline.
- Finish the runtime and consumer inventory, including lower-level wave consumers that do not mention `IMModes`.
- Record exact total-count, zero-count, ordering, explicit branch-selection and unavailable-sector behavior.
- Record representative F/G amplitudes, norms, signed metrics, reconstruction and source-projection results across numerical and analytical backends.
- Preserve current hydrostatic root-grid fixtures and manuscript plotted arrays/PDFs.

Acceptance: small executable target examples and the contract matrix identify every intended behavior change and every behavior that must remain invariant.

### Stage 1: establish the capable family foundation

- Add shared-background construction to the existing family problem path without duplicating physical state or validation.
- Put the manuscript wave convention and all associated metric/dual conversions in the wave family.
- Preserve capable analytical, spectral and finite-difference results, and provide direct special-family paths.
- Move descriptive normalization, metadata and no-solve selection into the capable bases as needed; compute norms lazily with explicit cache ownership and invalidation.

Acceptance: independent wave-equation, amplitude, normalization, derivative, endpoint, signed-metric and projection/reconstruction checks pass. A spectrum-only request does not evaluate profiles merely to choose a normalization. A normalization override or column selection does not repeat the eigensolve.

### Stage 2: consolidate requests and endpoint construction

- Apply the agreed total-count contract to scalar, analytical and bulk solves, including inertial pages.
- Preserve explicit requests for internal and subinertial inventories without automatic extra surface columns.
- Use explicit endpoint selection for canonical zero-APV responses and keep rotations and generalized-energy pairings available on their capable basis.
- Preserve empty inventories and exact auxiliary counts for neutral-label and critical-interval families.

Acceptance: the four-column APV example returns four columns without predicting the negative count; branch-specific, zero-count, insufficient-sector and repeated-page cases pass across supported backends.

### Stage 3: migrate integrations and remove the parallel layer

- Migrate hydrostatic root transforms and their returned basis, fixed-grid bridge, barotropic sector, root finding and qualification diagnostics.
- Migrate collections and repair F-derivative capability propagation.
- Relocate exact zero-depth wave construction and remove the analytical marker's dependence on the facade.
- Remove `IMModes` and the retired wrapper after all runtime dependants use the canonical family protocol. Avoid keeping a second permanent public construction route.

Acceptance: root-transform fixtures, fixed-grid recovery, retained-band qualification, neutral projections, endpoint rotations and collection tests retain their scientific coverage. Searches find no executable dependencies on the removed classes.

### Stage 4: migrate consumers and documentation

- Update manuscript scripts and any necessary WaveVortexModel consumers without changing unrelated dynamics or public behavior.
- Replace facade-specific tests with contract tests using the canonical API.
- Update README, canonical API help, migration documentation and CHANGELOG with actual before/after usage.
- Regenerate the seven manuscript figures once after the coherent code batch; compare numerical arrays and rendered figures.

Acceptance: manuscript scientific output is unchanged within declared numerical tolerances; any visible export differences are reported. Required WaveVortexModel construction, linear/source-projection, reconstruction, energy and restart checks pass or have a documented unchanged baseline limitation.

### Stage 5: final integration handoff

- Run affected focused suites, then the full provider `UnitTestsV2` suite because public wave behavior changes across construction paths.
- Run the relevant bounded consumer integration checks after the final provider changes.
- Run Code Analyzer once after the final coherent MATLAB batch, then rerun only files affected by repairs.
- Run whitespace, repository scope, local documentation/source-link, manifest and generated-artifact checks. There is no configured provider `docs:check` task; do not claim one ran.
- Maintain a verification ledger and report residual risks, actual missing assets, unchanged baseline failures and any explicitly deferred release qualification.

No source push, issue creation, PR creation or release publication is part of this local handoff. A future release must separately qualify the supported exported package graph.

## Tracking recommendation

Use this document as the local implementation checklist. The five design changes are coupled through the basis and metric contracts, so track one coordinated migration with sequential checkpoints rather than five independently executed refactors. Parallelize read-only review and disjoint family/consumer work after the shared contracts are fixed.

GitHub tracking is optional. If requested, one umbrella issue with the staged checklist is sufficient initially; create dependent implementation issues only where separate review or ownership helps. Posting requires current-task permission naming `JeffreyEarly/internal-modes` and the exact communication action, after confirming repository ownership under the shared workspace policy.

## Planning verification ledger

- Read the shared package-design, MATLAB-style, documentation-style and package-release guidance.
- Inspected provider working state, constructor/solver/basis contracts, direct facade/wrapper references, hydrostatic root-transform dependencies and analytical dispatch.
- Inspected manuscript callers and lower-level WaveVortexModel wave consumers.
- No MATLAB source or tests changed and no numerical gate rerun during planning.
- No external tracker communication performed.
