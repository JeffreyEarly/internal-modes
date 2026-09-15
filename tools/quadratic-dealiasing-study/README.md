# Shared quadratic-filtering score: issue #38

`assessQuadraticDealiasing` evaluates the agreed `none`, `fixedFraction`, and `effectiveBandwidth` policies from F/G values or coefficients in a caller-supplied common coordinate. The function returns per-mode masks; the consumer owns prefix selection, explicit counts, and production defaults. It performs no eigensolves or mode-pair products. See `Documentation/QuadraticDealiasing.md` for the mathematical contract.

Verified on 14 September 2026, Apple M4 Max, MATLAB R2026a Update 5. All 11 `IMQuadraticDealiasingTests` pass on a clean MATLAB path containing only this provider root and its `UnitTestsV2` folder. Source/test Code Analyzer and whitespace checks pass. GPT-6 Astra extra-high review found a full-coverage tiny-tail rounding issue; the correction and regression were reviewed again with no remaining findings. Full coverage inspects nonzero coefficients before squaring, and tail energies are summed directly.

The standalone `measureQuadraticDealiasing` script measures two dense synthetic polynomial channels, with one warmup and five repetitions of twenty calls, at `maxNumCompThreads(1)`. Raw measurements are in `provider-results.json`.

| Samples Q | Modes M | Median score time |
| ---: | ---: | ---: |
| 129 | 32 | 0.786 ms |
| 313 | 64 | 0.649 ms |
| 625 | 128 | 2.301 ms |
| 1249 | 256 | 9.191 ms |

These timings exclude evaluation of the caller's basis, linear assessment, prefix selection, and model construction. The smaller cases show startup/timing variation and should not be used to fit a complexity exponent. In particular, multiplying one isolated timing by a wavenumber count is not a measured constructor cost.

For values, the score costs O(M Q log Q); coefficient inputs cost O(M Q). Working storage is O(Q min(M,128)) plus O(M) reporting, excluding the caller's input arrays and FFT-library workspace. The Q=313, M=64 F/G inputs contain 320,512 bytes. This is an array-size calculation, not measured process RSS or a peak-memory bound. Dense basis evaluation, when required by a caller, must be counted separately and can be cubic if sample, polynomial, and mode counts scale together.

Reproduce with this provider on the MATLAB path, then add this study folder and run `measureQuadraticDealiasing("provider-results.json")`. Coordinate and sample-resolution adequacy for actual wave/APV modes, final defaults, and nonlinear behavior belong to the WVM integration and calibration issues #39–#40. The next target is #39: consume the score in existing common-grid preparation and verify constructor, persistence, and nonlinear-registration behavior.
