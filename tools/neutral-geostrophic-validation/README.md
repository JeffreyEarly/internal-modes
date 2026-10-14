# Neutral-stratification readiness validation

`runNeutralGeostrophicValidation` qualifies the observable balanced families before rerunning the coordinate-grid comparison. It does not run that study or change its archives.

From the repository root in MATLAB, reproduce the qualified joined-grid run:

```matlab
addpath(pwd,"tools/neutral-geostrophic-validation");
report = runNeutralGeostrophicValidation(discretization="piecewise");
```

For a focused diagnostic rerun, pass `caseNames="gpe-finite-1000"` and a separate `outputPath`; the report records the selected cases. The full readiness gate uses the default complete case list.

The joined-grid trial uses `discretization="piecewise"`. It splits at -400 and -200 m, keeping the total coefficient budget at 384 and 512 per solved field, with the same one-time 1024 retry. Points are allocated in proportion to interval depth, with at least 16 per interval. The report records the allocation and includes interface continuity in the equation residual. APV, generalized enstrophy, and zero APV enforce continuity of both F and G; MDA enforces G and G derivatives, with continuous F obtained by integration from the surface.

To retain both complete qualification runs:

```matlab
folder = "tools/neutral-geostrophic-validation";
globalReport = runNeutralGeostrophicValidation(discretization="global",outputPath=fullfile(folder,"global-readiness-report.json"));
piecewiseReport = runNeutralGeostrophicValidation(discretization="piecewise",outputPath=fullfile(folder,"piecewise-readiness-report.json"));
```

`neutralGeostrophicPiecewiseBasis` is a study-only generalized-eigenvalue solver. It uses the existing local spectral matrices and independently assembles the joined pencil; shooting roots are never used to find or refine its eigenpairs. Production factory and solver defaults are unchanged. The readiness flag applies to the discretization recorded in each report.

The current [handoff report](readiness-report.json) matches the [piecewise report](piecewise-readiness-report.json), with 12/12 cases qualified. The [corrected global report](global-readiness-report.json) retains the single-interval comparison; see [verification.md](verification.md) for the results table. Calling the runner without a `discretization` option retains its original global-grid default.

The default run checks the first 40 indexed modes, adding APV/MDA negative and constant sectors, and both zero-APV endpoint responses. It uses the prescribed neutral surface layer, manuscript and WVM endpoint choices, Earth rotation at 45 degrees, wavelengths 10 km and 1 km, and default, finite, and infinite generalized-enstrophy endpoint weights.

The independent reference uses the regular first-order equations. Positive roots come from endpoint-determinant shooting; the number of negative roots is checked against endpoint inertia. Negative modes use local transfer matrices and multiple shooting to preserve strongly localized structures. Neutral and constant-stratification intervals have exact propagators; variable intervals use `ode113` with tolerances `1e-11` and `1e-13`. Refined references reuse brackets found by the independent coarse scan, never production eigenvalues.

Physical-depth spectral references use 384 and 512 points, with one retry at 1024. Mode errors use a common positive physical-depth H1 norm for each stored field, including independent differentiation of the production representations. An exactly zero reference channel uses the companion-channel norm as its absolute scale; infinite equivalent depths are checked by exact zero inverse depth. Independent Clenshaw--Curtis quadrature uses 257 and 513 points on each profile interval, with one retry at 1025.

The MDA sampling check retains the full n+1-coefficient pressure antiderivative. It differentiates that stored polynomial independently of the differential equations; it also differentiates the stored G coefficients directly. Piecewise derivatives are evaluated from each stored local polynomial.

The study-only MDA solver uses Chebfun's `standardChop` at machine epsilon to remove unresolved local coefficient tails before differentiating twice. It pads the arrays back to the allocated pencil size and rechecks the original generalized-pencil residual. The report records cutoff degrees, coefficient perturbations, and the residual after cleanup. This representation cleanup changes neither the stratification nor the qualification tolerances. Mathematically exact constant sectors are represented exactly.

The JSON report records source SHA-256 identities, configuration, branch counts, mode labels, per-mode errors, refinement, quadrature convergence, equation/boundary residuals, and Gram checks. APV and generalized-enstrophy checks include both manuscript Gram identities; MDA checks its signed G Gram matrix. For zero APV, `maximumGramError` records the canonical endpoint-response identity error. JSON `null` equivalent depths are identified by `infiniteEquivalentDepth`; other nonfinite metrics indicate failure. A case remains `inconclusive` if any required reference, accuracy, or convergence gate fails. `readyForNeutralBalancedComparisons` is true only when every case qualifies and source identities remain fixed throughout the run.

This reference covers the observable balanced families. It does not construct independent neutral-label states with zero physical energy or establish completeness for those states.
