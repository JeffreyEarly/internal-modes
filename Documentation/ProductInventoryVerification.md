# Product inventory verification

This records #19 on `feature/sparse-product-inventories`. The implementation was developed against the #18 foundation, then integrated with exact bulk commit `ee7d340` before the final checks. The PR stack orders integration; #19's logical API dependency remains #18.

## Completed checks

- The final integrated batch passed all 48 tests: 17 inventory tests, one six-control frozen-fixture test, 11 shared projection tests, six collection tests, seven bulk solver tests, and six bounded-evaluation tests.
- The executable `ExamplesV2/BudgetedProductAssessment.m` passed. Its strict eight-column request rejects while retaining the request; the largest passing examined prefix remains a separate diagnostic.
- Code Analyzer reported zero findings for all three new classes and their methods, both test files, the example, and the documentation builder.
- Generated API documentation once after integrating the coherent source batch. Checked 929 generated Markdown files and all 863 InternalModes links, with zero missing targets.
- Verified the committed numerical fixture's SHA-256 and byte count against its readable manifest. The provider test runs without WVM classes or temporary study artifacts.
- Whitespace, package manifest/dependencies, changed-file scope, and absence of released snapshot changes checked before commit.

This repository has no configured `docs:check` task. Generator and explicit link checks are the available gates; no separate rendered-site or determinism verification is claimed.

## Scientific replay

The full pinned `cal-constant-17--fixed` WVM case reproduced 79,968 selected products, including 36,368 structural zeros. All 639,744 product-prefix comparisons agreed within 5.68e-14 after matching the original storage precision. The maximum aggregate prefix difference was 3.91e-9. See the [fixture provenance](../UnitTestsV2/Fixtures/README.md) and [manifest](../UnitTestsV2/Fixtures/manifest.json) for exact values, source versions, costs, reference discrepancies, and the source-version mismatch corrected during development.

Two additional [frozen replays](ProductInventoryReplays/README.md), exponential calibration and withheld pycnocline, passed all 639,744 and 3,779,072 product-prefix comparisons. Across the three cases, 5,058,560 comparisons passed. All structural-zero flags and nonzero counts matched; variable-stratification source channels and the original reference limitations were preserved.

The committed fixture provides six representative mixed controls; the full replay remains separately recorded evidence. Exact zero products, signed endpoint contributions, trigonometric same-family results, derivatives, coefficients, complex source pairings, and independent reference disagreement also have executable tests. Scientific source preparation was not timed in the initial constant replay, so it does not support an end-to-end speedup claim.

## Review corrections

- Reserve sparse ordinal anchors before expanding sign/column multiplicity; do not allocate a full factor Cartesian product just to reject a budget. A 100,000-column unrequested-tail test covers this distinction.
- Compare coverage by product and retained prefix, preserving a failure first exposed at an earlier untested prefix.
- Require the same scientific inventory, grids, shared prepared operators, and numerical condition guard before comparing coverage.
- Preserve a trustworthy measured failure even if another family has an inconclusive reference.
- Separate sampled projection failure from reference-system uncertainty.
- Require every requested retained ordinal and keep fixed boundary endpoints outside automatic modal truncation.
- Report available, examined, and omitted products at each requested count, including omissions within selected family rows.

These measurements do not resolve WVM's documented differentiated-pressure/vertical-momentum residual limitation, qualify arbitrary superpositions, adopt a default approximation, publish a provider release, or update WVM dependencies.
