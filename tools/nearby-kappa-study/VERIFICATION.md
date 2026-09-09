# Nearby-kappa research verification

- Completed the declared 96-case requested-pencil experiment and the separate 32-control product phase. Retained exact-repeat and common-metric triangle checks passed.
- Independent numerical review checked requested-κ equations, boundary signs, eigenvalue/frequency comparisons, positive subspace metrics, and fallback/oracle cost accounting. Added retained-band cutoff gaps so absent clustered-spectrum coverage remains explicit.
- Code Analyzer reported zero findings for `runNeighborKappaStudy.m`, `exportNeighborKappaResults.m`, and `assessNeighborKappaProducts.m` after bookkeeping cleanup. The manifest distinguishes the executed solve harness from its documented post-run preallocation/string/annotation cleanup; scientific operations and costs did not change.
- Integrated the tested provider commit `238f616`. Every recorded product-study provider/helper source hash and every recorded solve-study provider source hash matches the integrated checkout. No numerical rerun was needed for this unchanged-source integration.
- Checked result JSON parsing, local report links, whitespace, changed-file scope, unchanged package manifest/dependencies, and absence of released snapshot changes. This adds offline tools and reports only; no public class documentation was regenerated.

The study does not complete #10's accepted public solve-quality contract. Sharp-profile references, near-degenerate cutoff coverage, H1 quadrature qualification, independent source-field/physical-product qualification, and pressure/momentum residuals remain explicit limitations. The supplied evidence supports no default approximate reuse. Large continuous MAT artifacts remain in the recorded temporary output directories, with source and artifact hashes retained for review and reproduction.

The final whitespace check found CRLF endings in the profile CSV and trailing whitespace in the captured log. These text-format issues were normalized without changing numerical values; the check then passed.
