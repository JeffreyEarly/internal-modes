# Analytical Validation Summary

This note summarizes the V2 analytical-solution validation tests in
`UnitTestsV2`. The tests now exercise the canonical EVP stack:
`IMEigenvalueProblem`, `IMBoundaryCondition`, `IMInternalModes`,
`IMBasisSet`, and `IMInternalModesBasis`.

## Current Coverage

| Test suite | Analytical basis | Main purpose |
| --- | --- | --- |
| `IMConstantStratificationValidationTests.m` | `IMBasisSetConstantStratification` | Constant-stratification spectra, free-surface branch construction, hydrostatic `F` zero mode, and spectral coordinate smoke checks. |
| `IMExponentialStratificationValidationTests.m` | `IMBasisSetExponentialStratification` | Exponential Bessel-branch construction, free-surface branch construction, hydrostatic `F` null branch, and unsupported-boundary rejection. |

## Solver Validation

The numerical solver surface is `IMSolverSpectral` with
`coordinateKind="z"`, `"wkb"`, or `"density"`, plus
`IMSolverFiniteDifference`. The canonical tests compare assembly rows,
endpoint rows, mode selection diagnostics, and basis-set evaluation for these
solver paths.

## Unsupported Cases

| Area | Case | Expected result |
| --- | --- | --- |
| Constant analytical basis | unsupported endpoint condition | `IMBasisSetConstantStratification:UnsupportedBoundary` |
| Exponential analytical basis | nonzero surface location or unsupported endpoint condition | `IMBasisSetExponentialStratification:UnsupportedBoundary` or `IMBasisSetExponentialStratification:UnsupportedDomain` |
| Canonical mode selection | assessed positive metric and nonnegative form with raw negative discrete values | negative candidates are ignored during selection |
