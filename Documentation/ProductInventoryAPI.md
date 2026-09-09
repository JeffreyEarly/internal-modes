# Budgeted mixed-product assessment

`IMProductInventory` declares products without evaluating their bases. `fixedPlan` reserves a deterministic sparse selection and rejects insufficient budgets before any factor or output callback runs. `plan.assess` prepares the selected outputs, evaluates only selected factor columns on the required grids, and reuses those values and product pairings across retained counts. `IMProductAssessment.applyPolicy` makes a separate explicit decision and preserves the requested scientific counts. Its `requestedCount` refers to the inventory’s shared retained-count axis, which can map to two frequency-sign output columns per wave mode; it is distinct from the single-projection `requestedColumnCount` vocabulary. A trustworthy product failure rejects even when an unrelated family has inconclusive references; an unqualified numerical error alone cannot reject.

InternalModes owns numerical evaluation, signed projection, positive error norms, reference comparisons, and solver optimization. Callers own physical source channels, polarization/coefficient factors, valid Fourier interactions, endpoint injections, and model-state counts. The provider has no WVM dependency. The [basis/assessment contract](AssessmentAPI.md) remains the common foundation; new interfaces have no compatibility-wrapper requirement.

## Declare factors and products

```matlab
factor = IMProductInventory.collectionFactor(collection,id="page2-Gz",page=2,variable="G",derivativeOrder=1,countRole="retained",coefficient=@(z) 1+z);
inventory = IMProductInventory({factorA,factorB},products,{output});
plan = inventory.fixedPlan(interactionIds=["triad-a","triad-b"],prefixCounts=1:8,productBudget=10000);
result = plan.assess(grids,chunkSize=128);
decision = result.applyPolicy(quadraticTolerance=0.1,referenceTolerance=1e-4,requestedCount=8);
```

The executable [BudgetedProductAssessment.m](../ExamplesV2/BudgetedProductAssessment.m) supplies a complete trigonometric example.

A factor is a scalar struct with `id`, `family`, `labels`, `ordinals`, `frequencySigns`, `countRole`, `evaluate`, and nonempty `provenance`. The callback `evaluate(z,columns,referenceId)` returns an `nZ × nSelectedColumns` real or complex array. It must honor the explicitly selected columns, return finite values, and behave consistently throughout one execution. Evaluation has no hidden solve, fitting, or normalization change unless an explicit caller-owned preparation supplies and documents that operation.

`collectionFactor` binds an immutable collection snapshot, requested page, available variable, supported derivative, and a coefficient function. Its reference grids evaluate the same scientific modes. It does not turn a quadrature refinement into an independent solve. To assess independently solved modes, supply an explicit custom evaluator that selects and aligns a separately solved collection; preserve original mode labels, normalization, phase/orientation, and requested κ in provenance. Unsupported collection derivatives remain explicit errors; no derivative formula is invented from a family name.

`ordinals` identify positions along an explicitly varied retained-count axis. Multiple frequency signs may share one ordinal, with separate unique labels and `frequencySigns`. A `retained` factor must contain every ordinal from 1 through the maximum requested count, even when the requested counts are sparse. `fixed` factors keep every column at every count. Boundary endpoint factors must declare `columnKind="endpoint"` and `countRole="fixed"`; endpoints never become a modal prefix. APV, mean, inertial, and boundary counts can remain fixed independently of a varied wave count.

The product table has one row per supplied interaction/channel/ordered-factor/output combination:

| Column | Contract |
| --- | --- |
| `interactionId` | Caller-owned valid interaction identity; no Fourier geometry is inferred |
| `channel` | Physical or mathematical channel identity, supplied by the caller |
| `factorA`, `factorB` | Numeric indices into the factor cell row; product ordering is preserved |
| `output` | Numeric index into the output-recipe cell row |
| `coefficient` | Optional finite complex scalar, applied identically on every sample/reference grid |

Variable coefficient functions and derivatives belong in factor evaluators. Reuse the same factor identity for the same mathematical field/page/coefficient; different identities remain separate even if sampled values happen to agree.

## Selection and budgets

The initial policy is fixed sparse. For every count up through the requested band, retained factors select the low, middle, and cutoff ordinals `unique([1 2 ceil(count/2) max(1,count-2):count])`, restricted to the current count. The plan stores the cumulative union of pairs, preserving earlier tests. Fixed factors always retain all their columns, including signs and required endpoint coordinates.

The caller supplies already-valid interaction IDs; the provider preserves original inventory-row order for deterministic ties. The plan reserves every individual pair in every selected row, including products later found to be exactly zero. Budget rejection occurs entirely in metadata planning, before any factor evaluation, reference preparation, product projection, or independent reference solve.

`allProductsPlan` provides a bounded validation control over all pairs in the supplied interaction selection. It is useful for frozen-case comparisons and independent checks, not a claim of universal coverage. `result.compareCoverage(control,...)` requires the same scientific factor/output definitions, exact grids, the same numerical condition guard, and equal prepared projection/reference data for shared outputs. Callbacks must be deterministic; newly covered outputs use their identical declared recipe as provenance. It compares coverage per product **and retained count**, so testing a pair at a later count does not cover an earlier omitted count. It reports observed missed failures and false acceptance only when both references are qualified. Failures outside both evaluated selections remain unknown. No targeted/adaptive selection policy is implemented.

## Output and reference recipes

An output declares `id`, `family`, `labels`, nondecreasing `ordinals`, `countRole`, `prepare(grids)`, and explicit `provenance`. Preparation returns an `IMProjection` and a nonempty cell row of references. The projection must preserve every declared output label. A fixed output count is independent of the varied input count; retained output ordinals can group multiple frequency signs per count.

`grids` is a struct row of distinct `id` and finite real coordinate-column `z`. Its first ID is `sample`. Reference grids must have their own IDs. They may explicitly include endpoint coordinates alongside volume quadrature nodes; every reference operator and positive norm must use that same observation ordering. The provider neither guesses required endpoints nor drops supplied observations.

Each reference supplies:

| Field | Contract |
| --- | --- |
| `gridId` | A declared reference grid |
| `pairingMatrix` | Output-column-by-reference-observation signed physical pairing operator; complex values allowed |
| `normMatrix` | Real symmetric positive semidefinite product norm, including positive endpoint contributions |
| `targetGramMatrix` | Real symmetric signed coefficient system for this reference |
| `majorantGramMatrix` | Positive coefficient metric, validated on active output columns |
| `role` | Exactly one `primary`; optional `quadrature` and `independent` comparisons |
| `status` | `qualified`, `unverified`, or `inconclusive`; numerical agreement cannot override an unqualified status |
| `provenance` | Nonempty explicit reference construction, alignment, and convergence evidence |

The primary reference gives the reported product error. Other references compare retained coefficients and product normalization against the primary reference, using the assessed output's positive coefficient majorant. Quadrature and independent-solve discrepancies are recorded separately. An independent solve never silently replaces the assessed modes. Numerical comparisons do not manufacture external convergence evidence; `qualified` is an explicit caller assertion supported by its provenance.

Use `IMProjection.fromPairing` when a physical source dual is normalized by a supplied full-state coefficient system. WVM wave-source projection is one example: the energy denominator is not the sampled Gram of one source component. Such prescribed duals have no scalar synthesis or Gram/round-trip qualification. Separate basis-Gram and physical residual evidence remains necessary. Signed pairings are never replaced by absolute metrics; positive norms are never obtained by taking the absolute value of a signed quadratic form.

Exact zero references must also be zero on the sample grid. Nonzero products require strictly positive reference norms, without a small-denominator substitute. Inactive output coordinates require zero reference pairings. An explicit `minimumReciprocalCondition` guard can reproduce a caller's coefficient-system policy; use `1e-13` for the frozen WVM study. The default is zero, with singular systems still diagnosed. A failed sampled coefficient system gives a measured product rejection when references remain qualified. A failed reference system is inconclusive. Conditioning failures never silently change counts.

## Results and costs

`measurements` uses the common quantity/value/status vocabulary and records examined counts, structural zeros, limiting input labels, interaction, channel, and output. `evidence` retains every sampled product's input ordinals, labels, frequency signs, prefix errors, reference discrepancies, reference provenance, and condition guard result. `coverage` identifies selected and omitted inventory rows/interactions, per-family product counts, and explicit limits on interpretation. Its `prefixCoverage` table reports all available declared products, available products in selected family rows, examined products, and omissions within selected families and across the whole inventory at every retained count. These denominators are computed from metadata multiplicities without forming a dense Cartesian product. Exhaustive and superposition guarantees are false.

Execution caches are local to one call and keyed by the immutable inventory factor identity, selected columns, and exact named grid. No persistent global cache exists. New calls reevaluate their inputs. Factor calls and temporary product arrays obey `chunkSize`; cached selected factor arrays and returned product-error evidence occupy additional reported storage. Preparation, evaluation, and assessment times are separate. Memory counters report MATLAB `whos` array payloads and omit runtime, allocator, and opaque object overhead; they are not whole-process peak RSS.

The WVM frozen cases, signed/trigonometric controls, reference-convergence evidence, and costs originate in [study #400](https://github.com/JeffreyEarly/wave-vortex-model/issues/400), [PR #402](https://github.com/JeffreyEarly/wave-vortex-model/pull/402), and [advisory PR #405](https://github.com/JeffreyEarly/wave-vortex-model/pull/405). The documented differentiated-pressure/vertical-momentum residual limitation remains unresolved: product agreement is not complete model qualification. [InternalModes #10](https://github.com/JeffreyEarly/internal-modes/issues/10) retains ownership of solve-quality diagnostics and its existing milestone commitments. Provider release/export must precede any WVM dependency adoption.

The compact [frozen numerical controls](../UnitTestsV2/Fixtures/README.md) document their scientific provider version, full replay counts, reference discrepancy, measured execution costs, and source-version sensitivity. Their tests run without WVM.

Verification, frozen-case reproduction, review corrections, and validation limits are recorded in [ProductInventoryVerification.md](ProductInventoryVerification.md).
