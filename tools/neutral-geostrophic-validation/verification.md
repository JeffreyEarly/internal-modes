# Verification ledger

The task began with 134 modified or untracked entries from the rotation API and coordinate experiment. Those changes were preserved. Work stayed in the authoring checkout; the manifest, legacy v1 sources, WVM consumers, package snapshots, and coordinate-study archives were not changed.

## Solver tests

- Ran the v2 suite once: 438 passed and one high-resolution APV accuracy test failed.
- Corrected field and inverse-depth balancing, including row equilibration after physical field scaling. The high-resolution test passed after the change and again after final solver edits.
- Ran the final affected regressions: 135 passed, zero failed. This includes 17 paired-mode tests, six neutral supporting-family tests, existing APV/enstrophy/analytical tests, scalar endpoint/refactor tests, bulk wave tests, and wave-momentum tests.
- The merged full suite and final focused reruns contain 457 distinct passing tests, with no remaining failures. The final coherent rerun passed all 64 affected tests: 18 paired modes, eight staggered finite differences, 16 zero-APV tests, six neutral supporting-family tests, 12 orientation tests, and four independent-reference tests.
- The reference test correction uses one common F/G amplitude to remove the existing analytical family's separate quadrature normalization error; equivalent depths, physical shapes, and the independent reference Gram matrix retain their stated tolerances.
- The additional finite-difference smoke check exposed alternating modes in a collocated first-order discretization. Conservative staggered F/G assembly removes those branches. Tests verify second-order convergence against constant analytical solutions and neutral-cap spectral solutions, signed branches, the discrete Rayleigh identity, fully neutral and exact zero-inverse-depth families, and one zero-APV factorization per distinct wavenumber.
- After the localized-mode phase fix, all 30 paired-mode and orientation tests passed. Signed trapped modes whose surface traces are unresolved use a deterministic dominant-G phase, recorded in metadata; resolved surface traces retain the existing convention.

## Static and documentation checks

- Code Analyzer checked 33 changed MATLAB source/test/tool files. New solver, basis, and test files have no findings. A final targeted pass checked 13 files affected by the staggered implementation and final source freeze, with no production, test, or documentation-builder findings. Existing base-class advisory findings and two redundant cleanup-suppression notices in the validation runner remain; no scientific or functional analyzer issue was found.
- Generated API documentation, including both new classes and the changed factories and endpoint-response APIs. A local build-script workaround prevents the existing class-docs common-indent parser from truncating prose and mathematical markup. The successful rebuild was necessary after detecting this generator failure; the final phase and staggered-grid documentation was then regenerated as coherent follow-up batches after those changes.
- Restored 296 unrelated generated pages to their pre-task contents, including pre-existing local changes. Kept only task-related generated updates and the existing rotation documentation.
- Checked 910 generated local links. All new and task-affected links resolve. Two pre-existing tutorial links still point to absent `InternalModesSpectral/modesatfrequency` and `InternalModesSpectral/modesatwavenumber` pages; these were not changed.
- The documentation builder passed Code Analyzer with zero messages. Whitespace, tracked-file scope, manifest, and generated-artifact checks passed. This repository has no configured `docs:check` task.

## Neutral-layer readiness after the checker correction and split-grid trial

The MDA checker now retains all n+1 pressure-antiderivative coefficients. Paired modes and MDA G derivatives use their stored coefficients directly; no differential equation is substituted for a measured derivative. Four new sampling regressions cover the highest integrated coefficient, inconsistent F/G input, legacy integrated G, and piecewise derivative forwarding.

The study-only joined Chebyshev EVP splits at -400 and -200 m. APV, generalized enstrophy, and zero APV enforce F/G continuity. MDA enforces G/G-prime continuity and reconstructs continuous surface-referenced F by integrating the oversampled N2 G product. Its selected G series receives Chebfun standardChop cleanup at machine epsilon; eigenvalues are retained and the original pencil is rechecked. The maximum relative coefficient change is 5.60e-15, with a post-cleanup pencil residual below 3.89e-16. Actual stored-polynomial derivatives are used in every mathematical check.

Both full 12-case campaigns were rerun with identical original tolerances, independent shooting references, first 40 indexed modes plus required sectors, endpoint responses, and resolution budgets. Each report records 63 source identities, unchanged during its run and matching the saved files. Neither run raised a solver exception. See [README.md](README.md) for invocation and acceptance rules.

- [Corrected global report](global-readiness-report.json): 0/12 cases qualified after the 1024-point retry.
- [Piecewise report](piecewise-readiness-report.json): 12/12 cases qualified at 512 total coefficients per solved field; no 1024 retry was needed. The allocation is [307, 103, 102] across [-1000,-400], [-400,-200], and [-200,0] m, giving 510 unique physical nodes.
- [Latest handoff report](readiness-report.json): an identical copy of the qualified piecewise report, with `readyForNeutralBalancedComparisons=true` for that recorded discretization.

| Case | Global F H1 error (1024) | Piecewise F H1 error (512) | Piecewise equation/interface error |
|---|---:|---:|---:|
| apv-manuscript | 1.71e-08 | 6.26e-11 | 8.72e-09 |
| mda-manuscript | 1.64e-08 | 1.77e-11 | 2.26e-10 |
| apv-wvm | 1.37e-08 | 2.78e-11 | 5.41e-10 |
| mda-wvm | 1.37e-08 | 5.92e-11 | 5.94e-10 |
| gpe-default-10000 | 5.55e-08 | 3.22e-11 | 5.23e-10 |
| gpe-finite-10000 | 7.69e-08 | 2.24e-11 | 9.81e-10 |
| gpe-inactive-10000 | 1.38e-08 | 3.37e-11 | 3.79e-10 |
| zeroAPV-10000 | 1.05e-07 | 6.61e-14 | 1.26e-12 |
| gpe-default-1000 | 1.28e-06 | 1.30e-10 | 1.04e-09 |
| gpe-finite-1000 | 1.73e-06 | 2.04e-10 | 2.51e-09 |
| gpe-inactive-1000 | 1.38e-08 | 3.37e-11 | 3.79e-10 |
| zeroAPV-1000 | 2.26e-06 | 6.29e-12 | 3.75e-12 |

All piecewise cases pass independent-reference agreement, refinement, equivalent-depth and F/G H1 accuracy, differentiated equations and boundaries, interface continuity, quadrature convergence, and both applicable Gram identities. Worst-case errors are:

| Check | Maximum piecewise error | Required tolerance |
|---|---:|---:|
| Independent reference | 8.53e-10 | 1e-07 |
| Spectral refinement | 2.29e-10 | 1e-07 |
| Equivalent depth | 4.14e-11 | 1e-07 |
| F H1 | 2.04e-10 | 1e-07 |
| G H1 | 8.29e-11 | 1e-07 |
| Equations, boundaries, interfaces | 8.72e-09 | 1e-07 |
| Gram identities | 4.80e-10 | 2e-07 |
| Quadrature convergence | 3.12e-11 | 1e-09 |

All 14 final focused tests passed: four sampling, six piecewise-EVP, and four existing independent-reference methods. The new tests include constant analytical controls, exact generalized-enstrophy zero inverse depth, signed/constant sectors, canonical endpoint responses, neutral affine identities, and fixed point budgets. Code Analyzer checked all five affected MATLAB files: new helper and tests have no findings; the runner retains two pre-existing redundant cleanup-suppression notices. Whitespace and scope checks passed. Public solver/factory behavior, package metadata, generated website files, legacy v1, WVM consumers, package snapshots, and coordinate-grid archives were unchanged by this follow-up. No generated API sources changed, so no website rebuild or full production suite rerun was needed.

The results support the join-resolution explanation and provide qualified piecewise physical-z references for the subsequent neutral-layer balanced-mode comparisons. The trial also includes the documented MDA coefficient cleanup and oversampled pressure reconstruction. Wider solver adoption and the coordinate-grid comparison rerun remain separate tasks.

Two pre-existing generated documentation targets remain absent: `docs/classes/numerical-solvers/internalmodesspectral/modesatfrequency/index.md` and `docs/classes/numerical-solvers/internalmodesspectral/modesatwavenumber/index.md`, linked from `docs/tutorials/internal-modes-basics.md`. This follow-up did not modify those pages or links.
