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

The pages are organized around the objects that own each part of the calculation. Eigenvalue problems define the equations, solvers assemble and solve them, basis sets evaluate normalized modal variables, and analytical basis classes provide exact references for common stratifications.

## Choose a starting point

| Class group | Use it when |
| --- | --- |
| [`Core`](./core) | you want canonical eigenvalue-problem descriptors, endpoint conditions, scalar basis sets, or internal-mode basis sets |
| [`Solvers`](./solvers) | you want numerical solvers for physical, WKB-stretched, density-stretched, or finite-difference coordinates |
| [`Analytical bases`](./analytical-bases) | you want exact constant- or exponential-stratification basis sets |
| [`Supporting types`](./supporting-types) | you want shared normalization conventions |

## Shared notation

The reference pages use the standard vertical-mode variables:

| Symbol | Meaning | API name |
| --- | --- | --- |
| $$F_j(z)$$ | horizontal-velocity vertical structure function | `F` |
| $$G_j(z)$$ | vertical-velocity or density vertical structure function | `G` |
| $$h_j$$ | equivalent depth or eigendepth | `h` |
| $$N^2(z)$$ | buoyancy frequency squared | `N2` |
| $$f_0$$ | Coriolis parameter | `f0` |
| $$K$$ | horizontal wavenumber magnitude | `k` |
| $$\omega$$ | wave frequency | `omega` |
| $$z \in [-D, 0]$$ | physical vertical coordinate | `z`, `zDomain` |

## Reading the reference

- [`IMEigenvalueProblem`](./core/imeigenvalueproblem) defines the canonical scalar equation, endpoint conditions, normalization default, and grid-level diagnostics.
- [`IMInternalModes`](./core/iminternalmodes) translates standard `F` and `G` mode problems into canonical form.
- [`IMSolver`](./solvers/imsolver) is the abstract solver contract; concrete solvers such as [`IMSolverSpectral`](./solvers/imsolverspectral) and [`IMSolverFiniteDifference`](./solvers/imsolverfinitedifference) produce [`IMBasisSet`](./core/imbasisset) or [`IMInternalModesBasis`](./core/iminternalmodesbasis) objects.
- [`IMBoundaryCondition`](./core/imboundarycondition) stores scalar endpoint coefficients for the canonical boundary equation.
- [`IMBasisSetConstantStratification`](./analytical-bases/imbasissetconstantstratification) and [`IMBasisSetExponentialStratification`](./analytical-bases/imbasissetexponentialstratification) provide exact basis sets for common stratification profiles.
