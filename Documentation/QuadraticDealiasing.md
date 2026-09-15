# Quadratic dealiasing heuristics

`assessQuadraticDealiasing` measures individual mode shapes supplied by a caller and applies one of three inexpensive policies. The result is a heuristic filtering mask. It does not certify quadratic resolution, solve new modes, enumerate mode pairs, or select a model's retained count.

```matlab
% F and G have Q rows on a common, ascending WKB-Chebyshev grid.
% gridDegree describes the simulation grid; Q may be larger.
report = assessQuadraticDealiasing(F,G,gridDegree,quadraticDealiasing="effectiveBandwidth",energyFraction=.99,bandwidthFraction=2/3);
```

The caller must prepare both physical channels in one coordinate, such as the simulation grid's affine WKB coordinate. Value inputs are sampled at `xi=-cos(pi*(0:Q-1)'/(Q-1))`. Alternatively, `representation="coefficients"` accepts ascending Chebyshev degrees in that same coordinate. Coefficients taken directly from solvers with different native coordinates are not comparable.

For each nonzero channel, the effective degree is the smallest degree containing `energyFraction` of its Chebyshev-weighted energy, proportional to $2|a_0|^2+\sum_{n>0}|a_n|^2$. F and G are normalized independently, so changing either channel's nonzero amplitude leaves the criterion unchanged. A mode passes when the larger channel degree is at most `floor(bandwidthFraction*gridDegree)`. `[]` omits a channel, and exactly zero channels contribute degree zero.

`quadraticDealiasing="fixedFraction"` accepts the first `floor(retainedFraction*M)` columns, with an initial fraction of two thirds. Here M is the full linearly accepted inventory, before applying an explicit requested count. `quadraticDealiasing="none"` accepts every supplied column. Neither policy computes spectral diagnostics; its corresponding report fields contain NaNs.

The report includes an M-by-1 `accepted` mask, `effectiveDegreeF`, `effectiveDegreeG`, their maximum `effectiveDegree`, and `tailEnergyFractionF/G` above the permitted degree. It also records the policy, parameters, simulation degree, sample count, and mode count. The caller owns mode identities, ordering, contiguous-prefix selection, explicit-count handling, and any decision to require filtering for nonlinear evolution.

Transforming values costs $O(MQ\log Q)$; scoring supplied coefficients costs $O(MQ)$. Working storage is $O(Q\min(M,\mathtt{batchSize}))$ plus the $O(M)$ report, excluding caller-owned inputs. Batching uses base MATLAB's FFT and requires no parallel or optimization toolbox. Preparing the physical samples is a separate caller cost and can include dense basis evaluation. The fixed-fraction and no-filter paths skip the transforms and energy calculations, although input validation still scans supplied values.

`UnitTestsV2/IMQuadraticDealiasingTests.m` covers known spectra, weighted-energy thresholds, separate physical channels, large and small amplitudes, value/coefficient equivalence, native-coordinate conversion, bounded-batch equivalence, and empty or null inventories. The WVM calibration study selects production defaults; the defaults shown here are starting values rather than a nonlinear accuracy guarantee.
