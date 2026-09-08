# Unified assessment verification

This ledger records the #18 foundation on `feature/unified-mode-assessment`, based on `InternalModesEVP` at `9e2edb3`. It does not record completion of #19–#21 or qualification of WVM's physical model.

## Completed checks

- All 11 `IMProjectionTests` passed after prescribed-dual support, including complex pairings, signed normalization versus positive error metrics, inactive columns, rank deficiency, and prefix reconstruction.
- The final integrated batch passed `IMBasisAssessmentTests` (8), `IMDiscreteTransformPolicyTests` (10), `IMInternalModesDiscreteTransformPolicyTests` (8), `IMInternalModesDiscreteTransformTests` (8), `IMDiscreteTransformTests` (10), `IMBasisCollectionTests` (6), and `IMProjectionRecipeTests` (5).
- The affected certified-assessment and mean-density-anomaly suites also passed: `IMDiscreteTransformAssessmentTests` and `IMMeanDensityAnomalyModesTests`.
- `ExamplesV2/UnifiedBasisAssessment.m` executed successfully. It constructs scalar, aligned, APV, mean, wave, independent inertial, analytical, and both zero-APV boundary families.
- Code Analyzer reported zero findings for all new classes/tests and integrated MATLAB edits. The final parent pass covered the result/collection entry points, both existing transforms, both existing assessment routines, the new example, and the documentation builder.
- Generated API documentation once using `tools/build_website_documentation.m` and the local ClassDocumentation authoring package. Checked 897 generated Markdown files and 834 InternalModes links: no unresolved internal targets.
- Whitespace, source scope, unchanged package manifest/dependencies, and unchanged released snapshots checked before commit.

InternalModes has a documentation builder but no `buildfile.m` or configured `docs:check` task. The generator and explicit internal-link check above are the available checks; no separate determinism or rendered-site verification is claimed.

## Corrections found during verification

- Fixed scalar-prefix table orientation so single-band results have the same schema as multi-prefix results.
- Compared the old and new quadratic engines on the same explicit native reference quadrature, rather than confusing changes in reference integration with changes in the assessment kernel.
- Rejected inconsistent nonzero reference pairings on inactive columns instead of discarding them and accepting a false zero product.
- Defined product ordinals in selected assessment-column coordinates and tested reordered selections.
- Added the caller-prescribed physical dual because WVM's full-state energy normalization cannot be inferred from one scalar source component. Such a dual explicitly does not claim scalar Gram or round-trip assessment.

## Remaining series boundaries

Mixed-product inventory planning, budget preflight, and WVM frozen-case integration belong to #19. Shared bulk generation and bounded evaluation belong to #20. Nearby-κ reuse belongs to #21 and depends on relevant #10 solve-quality diagnostics. Independent reference evidence is always explicit. The WVM study's physical-residual limitation remains open, and provider release/export must precede a WVM dependency update.
