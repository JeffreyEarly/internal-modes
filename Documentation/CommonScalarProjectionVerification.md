# Common scalar adapter verification

This bounded #48 increment starts from `InternalModesEVP` at `e8fb1e8` (merged PR #65). It connects existing scalar physical recipes to the common selected/normalized basis and shared recovery engine. It does not complete #48's coupled state/source or augmented endpoint-observation adapters.

## Scientific controls and corrections

Independent author controls integrate APV F/G and MDA G forms on refined Gauss rules, including signed endpoint terms and positive majorants. Complex independent scalar inputs outside the retained span distinguish the physical dual from corrected recovery and least squares. Selection, reordering, normalization, empty inventories, prescribed DD APV, fixed collection pages, inactive G0 and deficient observations receive separate controls. Generic scalar tests retain a dense signed target, full positive majorant and string labels, without assuming a diagonal physical system.

Activity is explicit recipe metadata. The exact F-constant/G-zero sector is inactive in G for hydrostatic and APV Neumann/Neumann limits; the MDA constant-G sector remains active. No sampled magnitude or dimensional target threshold turns an active column into an inactive coordinate.

The independent mathematical audit found an inherited endpoint catalog defect: a relative activity threshold on canonical dimensional boundary coefficients could erase a genuine finite physical endpoint term. For example, free-surface APV with g=1e-16 and g0=.3 has a surface G coefficient nearly one, but the canonical coefficient comparison omitted it. APV G now uses the declared accelerations directly. MDA G likewise uses unit value functionals with coefficient g0/g or gd/g, avoiding reciprocal determinant and squared-functional intermediates. Generic boundary/solver classification is unchanged. Literal zero, inactive infinity and free-surface g0=-g reductions are explicit.

The adapter preserves existing source reference integration. Numerical forms use solver integration; analytical forms use a 1024-point trapezoidal rule. Independent Gauss comparisons expose the latter's finite error rather than treating analytic profile evaluation as exact integration. Those reference forms remain explicitly convergence-unverified. This increment does not claim a universal retained band or arbitrary-profile accuracy.

## Verification ledger

Verification uses MATLAB R2026b Update 1 on MACA64, one computational thread, with `matlab -batch` outside the local sandbox.

The first focused batch passed 77 of 83 methods. Six new controls failed because a hidden source recipe method was tested using public method lookup, causing the common adapter to choose its intrinsic empty/barotropic path. Source detection now uses the constructor's guaranteed source contract. The next run passed all existing lower APV/MDA regressions. Two new controls needed fixture/assertion corrections: a declared analytical background, and a dimensionless full-matrix error check in the positive majorant instead of an entrywise absolute bound on near-zero cross terms. Both failing methods then passed without changing the physical reference accuracy bound. The analytical G target/majorant operator error was 3.85e-6; numerical F/G target errors were at most 8.70e-16 in that control.

A separate numerical auditor passed **131 scientific/runtime controls**, including varying N², a different two-metre domain, signed accelerations, Gauss and irregular prescribed rules, the existing hydrostatic root rule, requested-empty and intrinsic barotropic inventories, a 1200-m dimensional normalization control, prescribed DD APV and failed recovery. Maximum measured relative errors were 9.72e-16 for H, 3.65e-16 for its majorant, 1.21e-15 for recovered round trips and 1.04e-15 for normalization covariance. The DD corrected round trip was 2.68e-16. Runtime instrumentation measured one shared engine preparation for a new scalar transform, none for binding prepared recovery, and no repeated EVP solve or weight fit. A tiny finite MDA coefficient of 1e-10 remained active and required its literal endpoint observation.

An independent mathematical auditor also passed four frozen-source controls. Selected signed APV physical projection agreed with refined Gauss integration to 2.43e-16, and corrected recovery agreed with independently assembled weak equations to 3.52e-16. Normalization covariance, MDA constant-G activity, APV zero-G inactivity and scaled endpoint coefficients passed. Independent framework review found no remaining required finding in the physical/recovery ownership, API, full-form transport or observation requirements.

Code Analyzer ran once across all **22 changed MATLAB files**. Four advisories in unchanged legacy `IMBasisSet` lines matched a separate base-revision check exactly; there were no new findings. Temporary API extraction succeeded for `IMModeBasis`, `IMScalarDiscreteTransform`, `IMProjection` and `IMBasisCollection`, producing **106 Markdown files** outside the repository. There is no configured `docs:check`; final source documentation checks supply that gate. No generated website content was edited.

The final full V2 integration passed **700 of 700 methods**, with no failures or incomplete methods. This includes all eight new common scalar methods and both new scalar binding methods. MATLAB production/tests were frozen and hashed before the full run; their hashes were verified afterward. Subsequent edits only clarified API comments for the common basis and prescribed-DD convenience properties; executable source lines were checked unchanged and the affected temporary API pages were regenerated.

Final source documentation checks covered four Markdown files, balanced code fences and six local links. Whitespace, repository scope, unchanged package metadata/dependencies and generated-artifact checks passed. The original unrelated manuscript worktree's tracked diff remained byte-for-byte unchanged. No released OceanKit snapshot was edited.

## Remaining work

Keep #48 open for coupled wave/state/source recovery, explicit derivative/companion endpoint observations, and other family adapters. MDA F and zero-APV coefficient-space forms do not acquire independent scalar physical projectors through this work. Cross-grid reference/retained-band qualification belongs to #58; external observation-layout fixtures to #60; provider grid calculus to #50; WVM adoption to #59. No release or package snapshot is produced.
