---
layout: default
title: Class documentation V2
nav_order: 5
has_children: true
permalink: /classes-v2
mathjax: true
---

# Class Documentation V2

This section is the API reference for building canonical eigenvalue problems, applying endpoint conditions, solving vertical modes, and evaluating basis sets in `internal-modes`.

The pages are organized around the objects that own each part of the calculation. Eigenvalue problems define the equations, solvers assemble and solve them, basis sets evaluate normalized modal variables, and analytical solution families provide exact references for common stratifications.

## Choose a starting point

| Class group | Use it when |
| --- | --- |
| [`Core`](./core) | you want canonical eigenvalue-problem descriptors, canonical or hydrostatic endpoint conditions, scalar basis sets, or internal-mode basis sets |
| [`Solvers`](./solvers) | you want numerical solvers for physical, WKB-stretched, density-stretched, or finite-difference coordinates |
| [`Analytical bases`](./analytical-bases) | you want exact constant- or exponential-stratification solution families, internal-mode bases, or SQG boundary modes |
| [`Supporting types`](./supporting-types) | you want internal-mode normalization conventions |

## Shared notation

The reference pages use the standard vertical-mode variables:

| Symbol | Meaning | API name |
| --- | --- | --- |
| $$F_j(z)$$ | horizontal-velocity vertical structure function | `F` |
| $$G_j(z)$$ | vertical-velocity or density vertical structure function | `G` |
| $$\lambda_j$$ | retained canonical eigenvalue | `eigenvalues` |
| $$h_j$$ | internal-mode equivalent depth or eigendepth | `h` |
| $$N^2(z)$$ | buoyancy frequency squared | `N2` |
| $$f_0$$ | Coriolis parameter | `f0` |
| $$K$$ | horizontal wavenumber magnitude | `k` |
| $$\omega$$ | wave frequency | `omega` |
| $$z \in [-D, 0]$$ | physical vertical coordinate | `z`, `zDomain` |

## Reading the reference

- [`IMEigenvalueProblem`](./core/imeigenvalueproblem) defines the canonical scalar equation, endpoint conditions, and grid-level diagnostics.
- [`IMInternalModes`](./core/iminternalmodes) translates standard `F` and `G` mode problems into canonical form and defines how eigenvalues map to equivalent depths.
- [`IMSolver`](./solvers/imsolver) is the abstract solver contract; concrete solvers such as [`IMSolverSpectral`](./solvers/imsolverspectral) and [`IMSolverFiniteDifference`](./solvers/imsolverfinitedifference) produce [`IMBasisSet`](./core/imbasisset) or [`IMInternalModesBasis`](./core/iminternalmodesbasis) objects.
- [`IMBoundaryCondition`](./core/imboundarycondition) stores scalar endpoint coefficients for the canonical boundary equation.
- [`IMHydrostaticBoundaryCondition`](./core/imhydrostaticboundarycondition) converts hydrostatic `F`/`G` endpoint laws into canonical boundary coefficients.
- [`IMConstantStratificationSolution`](./analytical-bases/imconstantstratificationsolution) and [`IMExponentialStratificationSolution`](./analytical-bases/imexponentialstratificationsolution) provide exact solution families for common stratification profiles.
