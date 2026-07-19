# Internal Modes V2 Discrete-Numerics Roadmap

This roadmap covers discrete-transform work built on numerical `IMBasisSet` objects. The mathematical reference is `vertical-modes-and-how-to-use-them/main.tex` in the OceanKit repositories workspace. V1 transform and projection classes remain separate and should be used only as temporary validation oracles.

## Available Now

| Capability | Current API | Result |
| --- | --- | --- |
| Scalar Galerkin transforms | `IMDiscreteTransform` | Stores `forwardMatrix` and `inverseMatrix`, applies `transformForward` and `transformBack`, and reports Gram, round-trip, and conditioning diagnostics. |
| Fixed-point weight fitting | `quadratureWeightsForPoints` | Fits nonnegative full-depth weights with the normalized Gram Frobenius objective, supports custom linear objectives, and compares fitted and geometric rules. |
| Mode-root point grids | `pointsFromModeRoots` | Returns both endpoints plus roots of selected column `nModes+1`, obtains an auxiliary mode when necessary, and uses solver-native spectral roots in physical, WKB, and density coordinates. |
| Fitted transform construction | `discreteTransform(z=...,nModes=...)` | Fits weights when they are omitted, returns the resulting transform, and optionally returns the associated `IMQuadratureWeightFit` provenance. |

The built-in weight objective stores only the independent upper-triangle Gram rows and applies the required `sqrt(2)` factor to off-diagonal rows, preserving the full normalized Gram Frobenius objective while reducing the least-squares system size.

### Current Scalar Workflow

```matlab
z = basisSet.pointsFromModeRoots(nModes=8);
[transform,weightFit] = basisSet.discreteTransform(z=z,nModes=8);
```

Callers can instead obtain or customize the fitted weights explicitly:

```matlab
[weights,weightFit] = basisSet.quadratureWeightsForPoints(z=z,nModes=8);
transform = basisSet.discreteTransform(z=z,weights=weights,nModes=8);
```

### Current Limitations

- `IMBasisSet.discreteTransform` transforms the solved scalar variable rather than a coupled internal-mode `F/G` family.
- Callers must still supply the physical sample points `z`.
- Scalar transform builders support endpoint inner-product terms only when they can be evaluated from sampled endpoint values; derivative traces are not inferred from arbitrary samples.
- `pointsFromModeRoots` currently requires solver-native spectral root support and is unavailable for finite-difference solvers.

## Next: Coupled Hydrostatic Transforms

Generalize the scalar machinery at `IMInternalModesBasis`, where `F`, `G`, and `h` represent one physical mode family.

- Use one shared physical point grid and weight vector with separate sampled matrices and metrics for `F` and `G`.
- Use the diagnostic `G` family to generate the canonical mode-root grid even when the EVP was solved in the `F` formulation.
- Handle the hydrostatic null mode explicitly: `F_0` participates in the `F` transform, while its aligned `G` column is bookkeeping only and is excluded from the `G` Gram system.
- Fit stacked normalized `F` and `G` Gram objectives with configurable block weights while preserving the custom-objective extension point and mode-pair provenance.
- Return separate forward transforms and quality diagnostics for the `F` and `G` components while retaining shared mode numbers, points, and weights.
- Start with rigid-surface and rigid-bottom hydrostatic modes. Generalized and active endpoint laws remain later work until their coupled discrete metrics are derived.

Acceptance requires the constant-stratification limit to reproduce DCT-I/DST-I behavior and exponential stratification to improve both component Gram errors relative to geometric weights.

## Then: Automatic Point Selection

Complete the one-call workflow by allowing the point grid to be omitted:

```matlab
[transform,weightFit] = basisSet.discreteTransform(nModes=8);
```

- Generate the mode-root grid, fit weights, and construct the transform internally.
- Preserve explicit paths for caller-supplied points and weights.
- Continue returning fit provenance so the convenient API exposes the selected points, constraints, optimizer result, fitted weights, and geometric comparison.
- Use the scalar workflow for `IMBasisSet` and the coupled workflow for `IMInternalModesBasis`.
- Never silently reduce the requested mode count. Report inadequate resolution or quality and let the user choose a smaller retained band.

## Later Milestones

### Retained-Band Diagnostics

- Compute normalized Gram errors, round-trip errors, condition numbers, weight ranges, and integral diagnostics as functions of retained mode count.
- For coupled modes, report `F` and `G` errors separately and provide a combined worst-component summary.
- Distinguish diagnostics for one fixed quadrature rule from diagnostics that refit weights at each retained count.
- Add an optional quality-policy helper that recommends the largest reliable prefix for user-supplied tolerances; do not make automatic truncation the default.
- Preserve named physical mode numbers in every diagnostic table and curve.

### Optional Mode-Root Point Optimization

- Initialize from `pointsFromModeRoots` and keep the physical endpoints fixed.
- Parameterize interior points so ordering and a configurable minimum spacing are maintained throughout optimization.
- Refit weights for each candidate grid and optimize retained-band Gram objectives rather than coefficient round-trip error alone.
- Support custom nonlinear objectives in addition to the built-in scalar and coupled Parseval objectives.
- Compare optimized points against mode-root points using identical weight constraints and report point displacement, conditioning, positivity, and both fitted and geometric-weight errors.
- Keep optimization optional; the mode-root grid remains the deterministic default.

### Generalized Metrics And Solver Support

- Add root hooks for finite-difference solvers only with refinement and convergence tests; do not treat unverified interpolation sign changes as spectral-quality roots.
- Keep analytical solution families separate from numerical basis sets; add analytical point-generation APIs to those families only when an exact construction is available and useful.
- Support derivative-dependent endpoint metrics through solver-provided trace matrices rather than estimating endpoint derivatives from arbitrary samples.
- Retain direct `IMDiscreteTransform` construction as the validated custom-matrix path. Add a basis-level custom metric builder only if repeated use cases justify a public API.
- Extend coupled transforms to generalized, free-surface, and active boundary conditions after their discrete endpoint terms and mode alignment are validated.

### Atlas And Observational Projection

- Store a finite atlas with normalization, quadrature rule, retained-band diagnostics, and provenance.
- Evaluate atlas fingerprints at imposed point or bin-averaged observation locations.
- Compare sampled observation Gram matrices with the appropriate partial-depth continuous Gram target.
- Select identifiable modal subsets using visibility, conditioning, and leakage diagnostics.
- Add noise and data-quality metrics, regularized projections, and uncertainty propagation as separate observational concerns.
- Couple to WaveVortexModel and add persistence only after the V2 transform contracts are stable.

## Development Rules

- Keep V2 tests and examples in `UnitTestsV2` and `ExamplesV2`.
- Validate constant and exponential stratification before adding broader profiles.
- Treat V1 comparisons as temporary dependencies and replace them with formulas or frozen self-contained fixtures before V1 removal.
- Keep scalar and coupled APIs mathematically explicit; convenience methods should compose the lower-level operations rather than create a parallel implementation.
- Add no bridge to `InternalModesBasis` unless a later integration milestone explicitly requires it.
