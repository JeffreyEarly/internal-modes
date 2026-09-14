# Wave eigensolver study

This study compares the production row- and column-scaled generalized solve
`eig(A,B)` with the algebraically equivalent inverse-pencil solve
`eig(decomposition(A)\B)`. The prototype always validates the returned modes
against the original generalized pencil and does not edit production behavior.

The frozen comparison revisions are WaveVortexModel
`c29c92805d2356cf0158b8b32834367156fa6dda` and InternalModes
`0e74261ca8c6ab25994fda2da0e9a495e2cd9228`. The live InternalModes study was
also run at the latter revision.

Run it in a fresh MATLAB process with the optimized package roots:

```matlab
config = fullfile(pwd,"wave-vortex-model-initialization-perf","tools","free-surface-initialization-study");
addpath(config)
configureInitializationStudy("round2-optimized");
study = fullfile(pwd,"internal-modes-initialization-perf","tools","wave-eigensolver-study");
addpath(study)
results = runWaveEigensolverCorrectnessStudy();
timing = measurePartialWaveEigensolver();
```

The timing columns are diagnostic only. Controlled alternating timing belongs
in the WaveVortexModel initialization study harness.

`runWaveEigensolverCorrectnessStudy` covers constant, exponential, and strongly
varying positive stratification; resolutions 69, 104, and 156; and zero, small,
intermediate, and large wavenumbers. It compares 64 retained labels,
equivalent depths, and sampled `F` and `G`, and evaluates residuals in the
original pencil. The aggregated output is in `dense-candidate-summary.csv`.

`measurePartialWaveEigensolver` uses the same three profiles and four
wavenumbers at resolutions 104 and 156. It warms both methods, alternates their
order for 15 repetitions, fixes the Arnoldi start vector, and sets tolerance
`1e-12`. `partial-timing-summary.csv` aggregates all 48 cases; each reported
mean is the mean of 12 per-case median ratios.

## Findings

Three dense standard-eigenproblem transformations were rejected. Factoring the
left matrix produced original-pencil residuals as large as `1.33e-8` and a
`1.51e-2` relative sampled-`G` difference. Exact bottom-constraint reduction
in coefficient space and physical-grid metric inversion were slower on median
and changed retained equivalent depths by as much as `3.16e-6` and `1.95e-6`,
respectively.

`measurePartialWaveEigensolver` measures a reduced generalized `eigs` solve.
In 15 alternating repetitions for three positive stratifications and four
wavenumbers, it was slower than full QZ at `nEVP=104`. At `nEVP=156`, the mean
partial/full ratios were `0.555` for 38 modes and `0.806` for 64 modes. All 48
partial solves returned flag zero. At `nEVP=156`, the largest original-pencil
residual was `7.8e-14`; maximum relative eigenvalue differences from full QZ
were `1.59e-10` for 38 modes and `3.21e-10` for 64 modes.

This is not yet a production recommendation. Convergence and Ritz residuals
validate returned eigenpairs but do not prove that a nonsymmetric collocation
solve did not omit a lower eigenvalue. A production experiment should restrict
the partial solve to high-resolution independent references, use deterministic
iteration, verify original-pencil residuals and contiguous modal root counts,
and fall back to full QZ on any failed check. The cost of that completeness
check must be included in an end-to-end benchmark.
Contiguous modal root counts are a plausible certificate under this
Sturm--Liouville problem, but their correctness and runtime cost were not tested
here.
