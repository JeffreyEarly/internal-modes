# Bulk wave API verification

This ledger records #20 on `feature/bulk-wave-bases`, stacked on the #18 foundation commit `1ac06d9`.

## Completed checks

- Seven new bulk solver tests, 43 existing EVP refactor tests, and 19 endpoint tests passed.
- Six new bounded-evaluation tests, 16 numerical zero-APV boundary tests, and 11 analytical tests passed.
- After integrating the #18 foundation, a final batch passed all 45 tests across `IMBasisAssessmentTests`, `IMBasisCollectionTests`, `IMBulkWaveModesTests`, `IMBulkBasisEvaluationTests`, `IMDiscreteTransformAssessmentTests`, and `IMInternalModesDiscreteTransformPolicyTests`.
- Code Analyzer found no issues in the new solver, evaluation, tests, or construction benchmark. The EVP class retains five pre-existing unused-self findings. The evaluation benchmark has two intentional temporary-result assignment findings in its timed loops; it verifies scalar agreement outside those loops.
- Construction benchmark ran five measured repetitions after warmup, separating repeated requests from independently deduplicated scalar solves. The evaluation benchmark ran three repetitions and found exact equality of full results; streaming checksum relative discrepancy was approximately 2.82e-13.
- Separate fresh MATLAB processes measured full and streamed evaluation maximum resident set size and macOS peak footprint. These include runtime and construction, and are explicitly separate from component-array sizes.
- Generated API documentation once after the coherent batch. Checked all 899 generated Markdown files and 836 InternalModes links, with zero missing targets. No configured `docs:check` task exists in this repository; no rendered-site or determinism check is claimed.
- Reviewed whitespace, changed-file scope, unchanged `resources/mpackage.json`, and unchanged dependencies and released package snapshots before commit.

## Review corrections

- Apply EVP diagnostic recovery before normalization, matching scalar basis evaluation even for custom recovery functions.
- Bound source-page evaluation for built-in numerical and analytical boundary bases; report custom callback allocation limits explicitly.
- Restrict shared solver preparation to exact built-in classes, so an unknown subclass with κ-dependent configuration cannot silently reuse incorrect setup.
- Preserve requested page positions in streaming callbacks, including duplicate pages.

The measurements and reproduction entry points are documented in [BulkWaveBases.md](BulkWaveBases.md). No nearby-κ approximation, model-state count policy, provider release, or WVM dependency change is included. The WVM study's physical-residual limitation is not resolved by these tests.
