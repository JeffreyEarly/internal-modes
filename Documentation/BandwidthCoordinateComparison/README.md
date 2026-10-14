# Trial bandwidth coordinate comparison

Implemented in the `internal-modes-evp` authoring checkout. The rotation API uses `rotationRate` and `latitude`; the trial `wkb` monitor is the square root of the local squared-frequency bandwidth. The physical EVPs retain their traditional rotation terms.

At Earth rotation and latitude 45 degrees, the certified prefixes are unchanged for the stratified profiles. This study does not demonstrate a resolved-mode gain. The neutral monitor is invertible and its EVP solves, but its first candidate mode fails the accuracy test at both candidate resolutions. Keep the change as a local trial; defer broader adoption until the neutral discretization and independent-reference accuracy are understood.

## Counts

Counts are consecutive certified prefixes from the first mode. A reference-limited prefix is a verified prefix, not a measured maximum. The machine-readable paired table preserves qualifications for the entire requested inventory. No count reached the candidate ceiling.

| Profile | Wavelength (km) | nEVP | Old N | New sqrt(B) | New minus old | Qualification |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| Constant | 10 | 64 | 29 | 29 | 0 | Assessed |
| Constant | 10 | 128 | 67 | 67 | 0 | Assessed |
| Constant | 1 | 64 | 29 | 29 | 0 | Assessed |
| Constant | 1 | 128 | 67 | 67 | 0 | Assessed |
| Standard exponential | 10 | 64 | 29 | 29 | 0 | Candidate cutoff verified; higher references inconclusive |
| Standard exponential | 10 | 128 | 45 | 45 | 0 | Reference limited; remaining modes inconclusive |
| Standard exponential | 1 | 64 | 29 | 29 | 0 | Candidate cutoff verified; higher references inconclusive |
| Standard exponential | 1 | 128 | 36 | 36 | 0 | Reference limited; remaining modes inconclusive |
| Weak deep exponential | 10 | 64 | 9 | 9 | 0 | Reference limited; remaining modes inconclusive |
| Weak deep exponential | 10 | 128 | 9 | 9 | 0 | Reference limited; remaining modes inconclusive |
| Weak deep exponential | 1 | 64 | 6 | 6 | 0 | Reference limited; remaining modes inconclusive |
| Weak deep exponential | 1 | 128 | 6 | 6 | 0 | Reference limited; remaining modes inconclusive |
| Neutral surface layer | hydrostatic G | 64 | unsupported | 0 | — | Inconclusive reference; first candidate fails; incomplete finite inventory |
| Neutral surface layer | hydrostatic G | 128 | unsupported | 0 | — | Inconclusive reference; first candidate fails; incomplete finite inventory |

For the neutral case, zero means no certified prefix. The first-mode F H1 errors against the 1024-point reference are 9.3975 at 64 points and 3.3999 at 128 points; relative equivalent-depth errors are approximately 4.6% and 14.6%. The 512/1024 reference comparison itself has first-mode F H1 error 1.030085e-7, narrowly above the required 1e-7. Neither a resolved-mode benefit nor accurate neutral low modes are established.

## Method and evidence

The domain is [-1000, 0] m, rotationRate is 7.2921e-5 rad/s, and latitude is 45 degrees north. Constant N=1e-2 s^-1 is the control; the exponential profiles are N=1e-2 exp(z/400) and N=3e-3 exp(z/300). The neutral profile is N²=0 above -200 m, has a quintic smoothstep to (5e-3)² between -200 and -400 m, and is constant below. The first three use fixed-wavenumber wave EVPs with surface Gz=G/h and bottom G=0; the neutral case uses hydrostatic G modes with both endpoints rigid.

Old and new candidates solve identical physical EVPs. The study-only LegacyWKBSolver forces only coordinate rotation metadata to zero. Independent physical-z references use 384/512 points, with one 512/1024 refinement for failing cases. Relative equivalent-depth and both F/G H1 errors must meet 1e-7 for references and 1e-5 for candidates. Up to nEVP-2 modes are assessed using one common 1025-point physical-z quadrature. Quadrature convergence was not separately studied.

- [Paired count table](paired-coordinate-comparison.csv) and [candidate qualifications](coordinate-comparison.csv).
- [Per-mode convergence errors](convergence-measurements.csv), [label matches](mode-matches.csv), and [reference validation](reference-validation.csv).
- [Common physical quadrature](physical-quadrature.csv), [run provenance](provenance.json), and [constant-grid control](constant-grid-control.json).
- Reproduce with `bandwidthCoordinateComparison(outputDirectory="<new-directory>")` after adding `ExamplesV2` and the package dependencies to the MATLAB path. Existing output directories are rejected.

## Verification ledger

- MATLAB R2026a on local Apple Silicon; batch execution outside the macOS sandbox.
- New rotation/coordinate tests: 9 passed. Full v2 suite after API migration: 414 passed, 0 failed, 0 incomplete.
- Zero rotation exactly reproduces the study legacy coordinate and passes the existing exponential map/derivative tests. Constant-profile normalized grids agree to roundoff at both candidate resolutions.
- Code Analyzer checked 76 changed MATLAB files against their baseline. No added findings remain; the study annotates intentional bounded result aggregation. Existing baseline findings are retained.
- API documentation rebuilt once with cached tutorials. Existing generated formatting and seven unrelated generator changes were restored to keep the diff scoped.
- All 116 local links from the 47 changed/new API pages resolve; all 829 generated site links resolve. Two pre-existing relative links in the unchanged legacy tutorial remain broken (`modesatfrequency` and `modesatwavenumber`). No configured docs:check task exists; no rendered-site or documentation determinism check is claimed. Full link and scope results are recorded in [verification.json](verification.json).
- No WVM consumers, legacy v1 source, package snapshots, package manifest, or publication was changed.

A [follow-up study of a Helmholtz-smoothed bandwidth monitor](../SmoothedBandwidthCoordinateComparison/README.md) uses the same physical cases and accuracy gates. It substantially improves the neutral first mode at 128 points, while the original reference-gated counts remain unchanged and weak-profile higher-mode diagnostics suggest a cost. The original results and provenance above are preserved.
