# Nearby-kappa feasibility: offline research only

This experiment supports [InternalModes #21](https://github.com/JeffreyEarly/internal-modes/issues/21). It changes no production API, acceptance default, or milestone. [Issue #10](https://github.com/JeffreyEarly/internal-modes/issues/10) still owns the eventual public solve-quality contract and remains outside the current V2 milestone.

The experiment tests a single-vector neighbor initialization for requested-pencil shift-invert Arnoldi, separately from unchanged neighboring-basis approximation. It is not a complete search of subspace acceleration methods. Continuous basis identity follows [#18](https://github.com/JeffreyEarly/internal-modes/issues/18), exact bulk construction uses [#20](https://github.com/JeffreyEarly/internal-modes/issues/20), and the bounded transform/product phase uses [#19](https://github.com/JeffreyEarly/internal-modes/issues/19).

Starting evidence includes [WVM #400](https://github.com/JeffreyEarly/wave-vortex-model/issues/400), [study PR #402](https://github.com/JeffreyEarly/wave-vortex-model/pull/402), and [advisory API PR #405](https://github.com/JeffreyEarly/wave-vortex-model/pull/405). Constant, exponential, and pycnocline profiles retain the study formulas. The sharper double-pycnocline profile is a declared new stress case, not frozen calibration evidence. The prior differentiated-pressure limitation is preserved; no result here constitutes complete model qualification.

## Reproduce

Add this authoring checkout and its study folder to the MATLAB path. Run each phase in a fresh process when measuring process peak memory.

```matlab
addpath(pwd);
addpath(fullfile(pwd,"tools","nearby-kappa-study"));
runNeighborKappaStudy("/private/tmp/im21-study/pilot",caseIds=[2 14 50 86],maxSeconds=120);
runNeighborKappaStudy("/private/tmp/im21-study/full",caseIds=1:96,maxSeconds=600);
exportNeighborKappaResults("/private/tmp/im21-study/full",fullfile(pwd,"tools","nearby-kappa-study","results"));
```

`maxSeconds` is a cooperative budget checked between cases; it does not interrupt baseline construction, a running solve, or artifact writes. Large MAT artifacts remain outside the repository. Compact configuration, case inventory, measurements, and provenance are retained alongside this README. The raw solver-phase summary preserves its original pending-product status; the later product phase has a separate summary and report.

## Declared experiment

The 96 cases combine four profiles, two surface conditions, three positive dimensionless kappa anchors, and four relative displacements, including an exact-repeat control. Bands 8/16 use 20 reference columns. Candidate solves use 64 Chebyshev coefficients; independent reference refinement uses 96/144, with at most one 216-coefficient escalation per case. WKB coordinates are fixed. Zero-wavenumber inertial modes are outside this neighboring-positive-kappa experiment.

Experimental reference guards are eigenvalue/frequency change of 1e-7, positive H1 subspace error of 1e-5, nondimensional off-grid G-equation residual of 1e-6, and normalized endpoint residual of 1e-8. Warm/cold Arnoldi uses tolerance 1e-10, at most 120 iterations, Krylov dimension at most 48, and 20 Ritz candidates. Additional candidate-to-full-64 guards use eigenvalue error of 1e-8 and H1 subspace error of 1e-6. These are research comparisons, not accepted solve-quality policies or universal accuracy guarantees.

The positive H1 comparison uses the integral of |G|²/D + D|G_z|² on a fixed 513-point trapezoidal grid. Both normalized/sign-aligned individual-mode errors and invariant-subspace errors are reported. Original candidate/reference labels remain intact. Numerical and analytical labels can differ; comparison uses eigenvalue ordering, with the labels retained as evidence. The reported spectrum includes internal gaps and the cutoff gaps between modes 8/9 and 16/17. Concentration and participation measurements describe the observed modes; a sharp profile alone does not establish localization or sensitive-subspace coverage.

References use independent spectral resolutions; constant/exponential cases also attempt analytical controls. ODE shooting is not included in this initial experiment. The H1 quadrature itself has not received an independent refinement qualification. Subspace agreement does not establish normalization, individual-mode rotation, transform, or product agreement.

## Retained state and costs

Each environment MAT file retains continuous numerical source, exact, reference, warm, and cold bases. Source and requested kappa are recorded separately; approximate sharing never masquerades as an exact `IMBasisCollection` mapping.

Costs distinguish requested construction, factorization, eigensolve, finalization, candidate validation, independent reference construction/validation, and measured full-solve fallback. The neighbor is assumed already solved; its acquisition cost is reported separately. The full-64 comparison oracle is constructed even when a warm candidate succeeds. The total including references charges this oracle rather than assuming a free online validator. Falling back to full-64 does not resolve an under-resolved reference.

The completed bounded transform/product phase is described in [product-assessment.md](results/product-assessment.md). It uses generic GG→G and explicitly L2 FG→F controls, and reports baseline, reuse-induced, and total error separately. It does not supply a complete WVM physical product inventory. No default approximation or WVM dependency update is authorized by this study; provider release/export precedes any consumer adoption.
