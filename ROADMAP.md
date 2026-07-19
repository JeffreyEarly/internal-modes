# Internal Modes V2 Discrete-Numerics Roadmap

This roadmap covers the discrete-transform work built on numerical `IMBasisSet` objects. The mathematical reference is `vertical-modes-and-how-to-use-them/main.tex` in the OceanKit repositories workspace. V1 transform and projection classes remain separate and should be used only as temporary validation oracles.

## Current Foundation

| Phase | Status | Result |
| --- | --- | --- |
| Scalar Galerkin transform | Complete | `IMDiscreteTransform` builds forward and inverse matrices from supplied points and weights and reports Gram, round-trip, and conditioning diagnostics. |
| Fixed-point weight fitting | Complete | `quadratureWeightsForPoints` fits nonnegative full-depth weights with normalized Gram objectives, supports custom linear objectives, and compares fitted and geometric rules. |
| Mode-root quadrature grids | Complete | `pointsFromModeRoots(nModes=...)` returns endpoints plus roots of selected column `nModes+1`, obtains an auxiliary mode automatically, and uses native spectral root finding for physical, WKB, and density coordinates. |

These phases intentionally establish scalar mechanics first. They do not yet provide the final coupled internal-mode transform or the one-call production workflow.

## Phase 4: Coupled Internal-Mode Transforms

Generalize the scalar machinery at `IMInternalModesBasis`, where `F`, `G`, and `h` represent one physical mode family.

- Build one shared grid and increment vector for both variables, with separate sampled matrices and metrics for `F` and `G`.
- Use the diagnostic `G` family to generate the canonical root grid even when the EVP was solved in the `F` formulation.
- Handle the hydrostatic null mode explicitly: `F_0` participates in the `F` transform while its aligned `G` column is bookkeeping only and is excluded from the `G` Gram system.
- Fit stacked normalized `F` and `G` Gram objectives with configurable block weights while preserving the current custom-objective extension point.
- Return separate forward transforms and diagnostics for horizontal kinetic energy and potential energy while retaining common mode numbers, points, and weights.
- Start with rigid-surface/rigid-bottom hydrostatic modes. Add generalized and active endpoint laws only after their coupled discrete metrics are derived.

Acceptance requires the constant-stratification limit to reproduce DCT-I/DST-I behavior and exponential stratification to improve both component Gram errors relative to geometric weights.

## Phase 5: Automatic Transform Construction

Make the common workflow require only a basis set and retained mode count.

```matlab
transform = basisSet.discreteTransform(nModes=8);
```

- Generate the mode-root grid, fit weights, and construct the transform internally.
- Preserve the existing explicit paths for caller-supplied points and weights.
- Preserve fit provenance and diagnostics so the convenient API does not hide the chosen nodes, constraints, optimizer status, or comparison with geometric weights.
- Use the scalar workflow for `IMBasisSet` and the coupled workflow for `IMInternalModesBasis`.
- Never silently reduce the requested mode count. Report inadequate resolution or quality and let the user choose a smaller retained band.

## Phase 6: Retained-Band Diagnostics

Turn transform assessment into a first-class API rather than a collection of individual scalar diagnostics.

- Compute normalized Gram errors, round-trip errors, condition numbers, increment ranges, and integral diagnostics as functions of retained mode count.
- For coupled modes, report `F` and `G` errors separately and provide a combined worst-component summary.
- Distinguish diagnostics for one fixed quadrature rule from diagnostics that refit weights at each retained count.
- Add an optional quality-policy helper that recommends the largest reliable prefix for user-supplied tolerances; do not make automatic truncation the default.
- Preserve named physical mode numbers in every diagnostic table and curve.

## Phase 7: Quadrature-Node Optimization

Optimize point locations only after mode-root grids and fitted weights provide a stable reference.

- Initialize from `pointsFromModeRoots` and keep the physical endpoints fixed.
- Parameterize interior points so ordering and a configurable minimum spacing are maintained throughout optimization.
- Refit weights for each candidate grid and optimize retained-band Gram objectives rather than coefficient round-trip error alone.
- Support custom nonlinear objectives in addition to the built-in scalar and coupled Parseval objectives.
- Compare optimized nodes against mode-root nodes using identical increment constraints and report node displacement, conditioning, positivity, and both fitted and geometric-weight errors.
- Keep optimization optional; the mode-root grid remains the deterministic default.

## Phase 8: Generalized Metrics and Solver Support

Extend the foundation only where a solver can provide numerically controlled traces and roots.

- Add root hooks for finite-difference solvers only with refinement and convergence tests; do not treat unverified interpolation sign changes as spectral-quality roots.
- Decide whether analytical basis objects should share the numerical basis interface or expose analytical roots separately.
- Support derivative-dependent endpoint metrics through solver-provided trace matrices rather than estimating endpoint derivatives from arbitrary samples.
- Add an expert custom metric-matrix path with explicit symmetry, rank, and target-Gram validation.
- Extend coupled transforms to generalized, free-surface, and active boundary conditions after their discrete endpoint terms and mode alignment are validated.

## Phase 9: Atlas and Observational Projection

Build partial-depth and noisy-observation workflows from a validated full-depth atlas rather than redefining the modes from each observation set.

- Store a finite atlas with normalization, quadrature rule, retained-band diagnostics, and provenance.
- Evaluate atlas fingerprints at imposed point or bin-averaged observation locations.
- Compare sampled observation Gram matrices with the appropriate partial-depth continuous Gram target.
- Select identifiable modal subsets using visibility, conditioning, and leakage diagnostics.
- Add noise and data-quality metrics, regularized projections, and uncertainty propagation as separate observational concerns.
- Couple to WaveVortexModel and add persistence only after the V2 transform contracts are stable.

## Development Rules

- Keep V2 tests and examples in `UnitTestsV2` and `ExamplesV2`.
- Validate constant and exponential stratification before adding broader profiles.
- Treat v1 comparisons as temporary dependencies and replace them with formulas or frozen self-contained fixtures before v1 removal.
- Keep scalar and coupled APIs mathematically explicit; convenience methods should compose the lower-level operations rather than create a parallel implementation.
- Add no bridge to `InternalModesBasis` unless a later integration phase explicitly requires it.
