# Mode creation API verification

Verified on 2026-10-15 against the authoring checkout. All local MATLAB commands used `matlab -batch` outside the sandbox. No release, source push, or tracker communication was performed.

## Focused checks

| Coverage | Passing distinct methods |
| --- | ---: |
| Shared background and generalized-energy configuration | 10 |
| Family creation, branch selection, automatic special sectors and bulk pages | 17 |
| Standard basis, normalization, derivatives and prescribed APV projection | 12 |
| Standard-basis collection integration | 6 |
| Neutral construction and numerical/analytical endpoint agreement | 6 |
| Neutral-label profiles, constant kernel and APV projection | 13 |
| Rotation-aware bandwidth, classical WKB and coordinate calculus | 13 |
| Existing bulk wave, count, MDA, neutral geostrophic, collection and projection regressions | 71 |
| **Total** | **148** |

Related edits were batched. Failing integration tests were corrected and rerun selectively. The final integrated creation/basis batch passed 40 methods; the subsequent critical-column edge test passed separately. Previously passing focused suites were not repeated without a relevant source change.

Independent wave-energy integration exposed aliasing from normalizing squared spectral profiles with the native collocation quadrature. The standard basis now uses oversampled Gauss quadrature in the numerical coordinate with the physical Jacobian. Independent energy checks pass for the surface, internal and subinertial sectors. Numerical and analytical mode identities, shared F/G scaling, signed generalized-energy signatures and observable neutral APV coefficients were checked separately.

## Static and scope checks

Code Analyzer covered all changed MATLAB files. New files and the modified neutral-label constructor have no findings. The existing exact-class advisory in `IMBasisCollection/evaluateChunks.m` and unused-argument advisory in `IMSolverSpectral/innerProductGrid` remain on unchanged statements.

Whitespace, local documentation links, Markdown fences, repository scope, package manifest and generated-artifact checks passed. The package manifest and declared dependencies are unchanged; `legpts` uses the already declared Chebfun dependency. Released OceanKit snapshots and generated website files were not changed. This repository has no configured `docs:check` task; local source/link checks were used. Website publication is outside this task.

## Limits of this verification

No full-suite, clean exported-package, production WaveVortexModel or shared-grid study qualification is claimed. The standard basis supplies the neutral-label prescribed APV dual; generic wave projection and transform construction remain available through the existing lower-level APIs. Analytical backends report unsupported families and frequencies explicitly. Declared exact neutral/critical intervals are checked for sampled contradictions, but sampling cannot prove an arbitrary function is constant throughout an interval.
