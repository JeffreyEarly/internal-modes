# Basis normalization verification

This increment implements #47 from `InternalModesEVP` at `72d754e` (merged PR #63). It adds an immutable common-basis normalization operation and repairs selected-column normalization integration. It does not add a recovery policy, observation adapter, general rotation API or WVM dependency.

Verified on 2026-11-02 using MATLAB R2026b Update 1, MACA64. Every MATLAB command used `matlab -batch` outside the local sandbox with one computational thread. Applicable package, MATLAB, caching, documentation and focused agent guides were read.

## Numerical correction

The former intrinsic/analytical integration rule chose its order from the current retained-column count. Selecting analytical harmonics 160 and 137 before requesting G depth normalization reduced the Gauss rule to 128 points. Relative scale errors were approximately -2.7% and +0.11%, compared with normalization before selection.

The basis now captures its original source-inventory count as immutable implementation state. Selection preserves that count and therefore the normalization integration rule. Analytical and interval-supported high-harmonic regressions check exact amplitudes, selection commutation and dropping the highest source column. No quadrature-weight fit or new physical EVP solve is introduced. This preserves the rule across selections; it does not prove arbitrary input functions are resolved by that rule.

## Author controls and independent audit

The basis-operation author supplied 12 passing test methods covering fresh construction across all family rules, signatures, null/nonfinite errors, profiles and derivatives, selected/reordered/empty inventories, immutable records, same-name reuse, and absence of physical solve/weight-fit entry points. A numerical covariance author supplied 7 passing methods in a separate test class for physical and recovery operators. Production changes and test expectations were independently reviewed.

The independent read-only runtime audit passed **357 distinct assertion controls**. Extra references included APV, inertial and hydrostatic-G families, exponential analytical and finite-difference sources, signed MDA energy, and an adversarial intrinsic fixture whose first column contains harmonic 173. Dense complex coefficient scales were checked under least-squares, corrected-dual and extended-span policies. Failed recovery retained its separately available physical dual; empty roots and generalized-enstrophy requests retained their contracts.

Covariance checks use actual wave energy and signed MDA pairings, off-span complex inputs, fixed physical-depth and existing root F/G adapters, and the prescribed DD APV dual/majorant. They preserve observation and fitting factors, scientific activity and coefficient-scale meaning. The neutral duality target remains identity; physical coefficient forms transform by congruence.

Root covariance requires the original recovery policy and extended span. The default truncated-collocation control transports retained coordinates while preserving omitted latent coordinates and the public prepared full-span specification; fixed points and weights alone would not preserve that recovery. No new transform-normalization API is implied.

## Final integration and static gates

- The final full V2 run passed **666 of 666 methods**, with no failures or incomplete methods. This includes all 19 new methods and the two previously corrected baseline test fixtures. Production and tests were frozen before this run.
- Code Analyzer ran once across all four changed MATLAB files and reported zero findings.
- Temporary API documentation extraction ran once for `IMModeBasis`, producing 31 Markdown files including the new `withNormalization` entry. No generated website files were changed.
- Source Markdown fences and 24 local links passed checks; existing legacy README fences were preserved. There is no configured `docs:check` task in this repository, so source checks and temporary API extraction supply the relevant documentation gates.
- Whitespace, repository scope, unchanged package manifest/dependencies and generated-artifact checks passed. The unrelated original dirty manuscript worktree's tracked diff remains unchanged. No released OceanKit snapshot was edited.

## Scope of the evidence

Normalization divides all related common-family profiles and derivatives by one positive per-column magnitude; physical signed pairings remain separate from the current normalization signatures. Labels, spectra, horizontal roots, original construction requests/evidence and selection maps remain unchanged. Only the resolved normalization record changes. Unsupported or null rules fail without removing columns.

The tests establish supported coordinate covariance in representable ranges with consistent policy, spans, activity and scales. They do not establish arbitrary-function integration convergence, a larger retained band, universal grid recoverability or new family/source projection semantics. Those qualification boundaries remain explicit in the [normalization guide](BasisNormalization.md).
