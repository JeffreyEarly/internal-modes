# Initial nearby-kappa feasibility results

The bounded experiment found **no demonstrated advantage from the neighboring initial vector**. Warm and deterministic cold shift-invert Arnoldi both used a median of 49 operator applications. Both were faster than the full dense solve in raw solver time, but validation and reference costs exceeded that saving. These findings support retaining exact bulk solves while the solve-quality and product-assessment work proceeds. They do not establish that every possible subspace acceleration method is unhelpful.

All 96 cases completed against provider commit `ee7d340`. This is partial evidence for [#21](https://github.com/JeffreyEarly/internal-modes/issues/21): complete physical product inventories, independent physical pressure residuals, and an accepted [#10](https://github.com/JeffreyEarly/internal-modes/issues/10) diagnostic contract remain outstanding. No production approximation or public acceptance policy was added. The full experiment and provenance are described in [README.md](README.md).

## Accurate requested-pencil solves

Every warm candidate reproduced the first 16 columns of the requested full-64 discretization within the declared experimental comparison guards. The largest relative eigenvalue difference was 2.70e-11; the largest positive H1 subspace error was 3.54e-10. This establishes agreement with that discrete solve, not accuracy of its approximation to the continuous problem.

Only 32 of 96 independently refined references met all experimental reference guards. Only 24 cases also passed the candidate/reference guards. The other 72 cases require the full-solve fallback and retain unresolved-reference or candidate-resolution status; the fallback does not repair inadequate spectral resolution.

| Profile | Surface | References meeting experimental guards | Cases requiring fallback |
|---|---|---:|---:|
| Constant | Rigid | 12/12 | 0/12 |
| Constant | Free | 8/12 | 4/12 |
| Exponential | Rigid | 12/12 | 8/12 |
| Exponential | Free | 0/12 | 12/12 |
| Pycnocline | Rigid | 0/12 | 12/12 |
| Pycnocline | Free | 0/12 | 12/12 |
| Double sharp pycnocline | Rigid | 0/12 | 12/12 |
| Double sharp pycnocline | Free | 0/12 | 12/12 |

The low-kappa constant/free pilot control illustrates why the measurements remain separate. At the 216-coefficient reference resolution, its eigenvalues agreed with the analytical control to 4.18e-10 and its H1 subspace to 7.66e-11. The normalized endpoint residual was 5.65e-10 and the matrix backward residual was 6.98e-16. Nevertheless, the independently sampled, differentiated G-equation residual reached 1.11e-5, exceeding the declared 1e-6 guard. No threshold was relaxed. This is consistent with the study's differentiated external-mode sensitivity and prevents presenting excellent eigenvalue/shape agreement as complete physical qualification.

The constant/free analytical external mode carries label -1 while the numerical family carries label 1 for the corresponding eigenvalue. Both original label sets remain in the MAT records; comparison follows eigenvalue ordering. No label was silently rewritten.

## Unchanged neighboring bases

An unchanged neighboring basis can have a nearly identical span while its eigenvalues, frequency factors, or F scaling differ. Constant rigid G spans were essentially kappa invariant in this experiment. That observation does not justify sharing the entire wave basis or transform.

The table reports the maximum **total difference from the independently refined requested-kappa reference** among the eight reference-qualified cases at each offset, for the 16-column band. It includes the original full-64 discretization error; it is not a pure reuse-error estimate. Exact-repeat controls show that baseline floor. The completed [bounded operator assessment](results/product-assessment.md) separately measures exact baseline error, the incremental neighbor-minus-exact effect, and total referenced error.

| Relative kappa offset | Eigenvalue difference | H1 subspace difference | F difference |
|---:|---:|---:|---:|
| 0 | 2.16e-8 | 4.07e-9 | 1.56e-8 |
| 0.001 | 1.82e-3 | 2.90e-4 | 1.82e-3 |
| 0.01 | 1.80e-2 | 2.91e-3 | 1.83e-2 |
| 0.1 | 1.60e-1 | 3.04e-2 | 1.91e-1 |

Across the same qualified cases, the largest frequency difference was 9.14e-2. Proximity in kappa and small Gram or subspace errors therefore cannot be used as interchangeable checks of physical accuracy.

## Coverage and limitations

The free-surface external modes show measured spatial concentration: an upper-quarter H1 fraction above 0.995 and a participation fraction as low as 0.093. These are measurements of the returned reference modes, not labels inferred from the profile name.

The matrix did **not** encounter near-degenerate retained-band cutoffs. The smallest relative gap at either the 8/9 or 16/17 cutoff was 0.102; the smallest reported internal gap was about 0.048. Thus this study covers an external-mode transition and concentrated modes, but does not qualify reuse around tightly clustered eigenspaces. The sharp-profile cases also lacked qualified references within the bounded refinement sequence. Those gaps must remain explicit in any follow-up decision.

The H1 comparison uses a fixed 513-point positive quadrature and has not independently qualified that quadrature through refinement. Individual-mode errors and subspace errors are both retained, but there is no general cluster-tracking or branch-completeness guarantee. Exact analytical controls were attempted for constant/exponential profiles; an independent ODE shooting implementation was not added here.

## Costs and retained evidence

These are medians of per-case measurements, so their components should not be added as if they came from one representative case.

| Measurement | Median |
|---|---:|
| Full requested solve, including amortized shared preparation | 3.696 ms |
| Warm shift-invert factorization/solve/finalization | 2.317 ms |
| Cold shift-invert factorization/solve/finalization | 2.235 ms |
| Candidate, stored-oracle validation, and required fallback | 25.843 ms |
| Candidate plus all reference/full-64 oracle acquisition and validation | 185.352 ms |
| Neighbor acquisition, reported separately | 5.484 ms |

The 25.843 ms quantity assumes already stored reference/full-64 oracles for successful candidates. The 185.352 ms quantity charges independent reference construction, reference validation, baseline validation, and the full-64 comparison oracle even when a candidate succeeds. Both assume an already solved neighbor; its acquisition cost is reported separately. The full solve was retained for every case to make the experiment reviewable. Neither accounting convention demonstrates an online validation method that is free of those reference costs.

The research computation took 52.62 s; the fresh MATLAB process took 55.54 s. Maximum resident set size was 967,720,960 bytes and peak memory footprint was 602,703,696 bytes. Large MAT artifacts total 404,729,213 bytes and remain under `/private/tmp/im21-study/full`; they preserve continuous source, exact, reference, warm, and cold bases. The repository contains compact CSV/JSON evidence and hashes in [artifact-manifest.json](results/artifact-manifest.json), along with [process-resources.log](results/process-resources.log).

The full run's raw log includes Code Analyzer notices subsequently addressed by case-inventory preallocation, string-comparison syntax, and annotation cleanup. The numerical measurements and cost formulas were not changed after execution. Both executed and current harness hashes are recorded in the manifest.

## Transform/product follow-through

The #19 inventory API completed 32 generic controls with 10,464 reserved output products. Fifteen met the additional experimental reference guard; seventeen remained inconclusive. At a 10% κ offset, constant rigid-lid G reuse changed coefficients by about 4.9e-12, while unchanged F in the declared L2 control changed coefficients by about 0.146. Exponential rigid-lid G reuse also produced about 0.146 error. Exact-repeat reuse differences were zero and per-product common-metric triangle checks passed. See the [full product results](results/product-assessment.md) for pairing definitions, coverage, source-field limitations, and compact data.

This supports retaining exact bulk generation as the current path. A future approximation would need an accepted solve-quality contract, references that resolve the sharp cases, tightly clustered spectral coverage, and the relevant caller-supplied physical products. The investigation has not authorized default approximate sharing.
