# Wave eigensolver study

This study compares the production row- and column-scaled generalized solve `eig(A,B)` with three dense standard-eigenproblem transformations and a reduced generalized `eigs` solve. Every timed path requests eigenvectors and eigenvalues, reconstructs the vectors in the solver's native representation, and reports preparation, solve, reconstruction, and total costs separately. The prototypes validate returned modes against the original generalized pencil and do not edit production behavior.

The archived comparison used WaveVortexModel `c29c92805d2356cf0158b8b32834367156fa6dda` and InternalModes `0e74261ca8c6ab25994fda2da0e9a495e2cd9228`. Its timing ratios are superseded because the full and partial paths did not use matched output contracts and excluded required preparation and native-vector reconstruction costs. The old files remain as `dense-candidate-summary-superseded-output-contract.csv` and `partial-timing-summary-superseded-output-contract.csv`.

Run the corrected studies in a fresh MATLAB process with the optimized package roots:

```matlab
config = fullfile(pwd,"wave-vortex-model-initialization-perf","tools","free-surface-initialization-study");
addpath(config)
configureInitializationStudy("round4-full");
study = fullfile(pwd,"internal-modes-initialization-perf","tools","wave-eigensolver-study");
addpath(study)
[denseSummary,partialSummary,correctnessResults,partialResults] = refreshWaveEigensolverStudySummaries();
```

`refreshWaveEigensolverStudySummaries` runs both corrected studies and writes `dense-candidate-summary.csv` and `partial-timing-summary.csv` only after both result tables have been assembled. The checked-in canonical summaries were refreshed with the active `round4-full` package-root configuration on 2026-09-13.

`runWaveEigensolverCorrectnessStudy` covers constant, exponential, and strongly varying positive stratification; resolutions 69, 104, and 156; and zero, small, intermediate, and large wavenumbers. It compares 64 retained labels, equivalent depths, and sampled `F` and `G`, and evaluates residuals in the original pencil. For every dense candidate, its total and the generalized baseline total both begin with assembled `A,B`, request `[V,D]`, and end after native-vector reconstruction.

`measurePartialWaveEigensolver` uses the same three profiles and four wavenumbers at resolutions 104 and 156. It warms both complete output paths, alternates their order for 15 repetitions, fixes the Arnoldi start vector, and sets tolerance `1e-12`. The total full timing includes scaling, `[V,D] = eig(...)`, and native-vector reconstruction. The total partial timing additionally includes exact bottom-constraint reduction and reduced-pencil formation before `[V,D,flag] = eigs(...)`.

`runWaveEigensolverPartialStudy` is a smaller exploratory run for one exponential profile. It now uses the same full production baseline, complete partial path, phase boundaries, correctness checks, and explicit failure reporting.

The dense-candidate and exploratory timings remain single-pass diagnostics, and the generalized baseline always runs before each dense candidate. Only `measurePartialWaveEigensolver` provides warmed, alternating, repeated timings suitable for comparing the measured full and partial paths.

## Findings

The archived correctness checks rejected the three dense standard-eigenproblem transformations. Factoring the left matrix produced original-pencil residuals as large as `1.33e-8` and a `1.51e-2` relative sampled-`G` difference. Exact bottom-constraint reduction in coefficient space and physical-grid metric inversion changed retained equivalent depths by as much as `3.16e-6` and `1.95e-6`, respectively. Their archived speed ratios are superseded.

The archived partial/full ratios are not comparable to production because they timed eigenvalue-only calls and gave the partial method its reduction and reconstruction work for free. In the corrected controlled timing, the mean partial-total/full-total ratios are `1.110` and `1.166` for 38 and 64 requested modes at `nEVP=104`, and `0.434` and `0.615` at `nEVP=156`. These totals include all measured scaling, reduction, solve, and native-reconstruction work.

All 48 corrected partial cases returned flag zero. The maximum relative eigenvalue difference from full QZ was `3.21e-10`, and the largest original-pencil residual was `8.68e-13`. The dense single-pass diagnostic ratios were `3.23` for `gridMetric`, `1.93` for `inverseA`, and `1.85` for `reducedMetric`; those figures are context rather than controlled performance estimates. The dense correctness results remain unchanged in substance.

This is not yet a production recommendation. Convergence and Ritz residuals validate returned eigenpairs but do not prove that a nonsymmetric collocation solve did not omit a lower eigenvalue. A production experiment should restrict the partial solve to high-resolution independent references, use deterministic iteration, verify original-pencil residuals and contiguous modal root counts, and fall back to full QZ on any failed check. The cost of that completeness check must be included in an end-to-end benchmark. Contiguous modal root counts are a plausible certificate under this Sturm--Liouville problem, but their correctness and runtime cost were not tested here.
