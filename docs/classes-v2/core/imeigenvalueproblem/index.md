---
layout: default
title: IMEigenvalueProblem
has_children: false
has_toc: false
mathjax: true
parent: Core
grand_parent: Class documentation V2
nav_order: 1
---

#  IMEigenvalueProblem

Describe a vertical-mode generalized eigenvalue problem.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMEigenvalueProblem</code></pre></div></div>

## Overview

`IMEigenvalueProblem` is the solver-independent contract for a
vertical-mode eigenvalue problem. A solver owns the stratification,
physical domain, coordinate mapping, and derivative matrices. The EVP
owns the physical constants, differential operators, boundary laws,
inner-product weights, normalization rules, equivalent-depth mapping,
and mode-index policy.

Assembly combines those responsibilities in the solver native basis:
$$Aq_j=\lambda_jBq_j,\qquad h_j=\mathrm{hFromEigenvalue}(\lambda_j).$$
The retained columns become an `IMBasisSet`. The basis set evaluates the
solved variable and its linked diagnostic variable; `F` and `G` are not
independent mode families. In a `G` formulation,
$$F_j=h_j\partial_zG_j,$$
and in an `F` formulation,
$$G_j=-gN^{-2}\partial_zF_j.$$

Use the static factories for standard wave and hydrostatic mode
problems. Use the constructor directly when defining a custom operator,
boundary, inner-product, or normalization contract.

```matlab
N2 = @(z) 1e-5*ones(size(z));
solver = IMSolverSpectral(N2=N2, zDomain=[-1000 0], nEVP=64);
evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=1e-4);
basisSet = solver.solveEVP(evp, nModes=4);
z = linspace(-1000, 0, 128).';
G = basisSet.G(z);
h = basisSet.h;
```




## Topics
+ Create standard EVPs
  + [`hydrostaticFModes`](/internal-modes/classes-v2/core/imeigenvalueproblem/hydrostaticfmodes.html) Create the geostrophic hydrostatic `F`-mode EVP.
  + [`hydrostaticGModes`](/internal-modes/classes-v2/core/imeigenvalueproblem/hydrostaticgmodes.html) Create the hydrostatic `G`-mode EVP.
  + [`waveModesAtFrequency`](/internal-modes/classes-v2/core/imeigenvalueproblem/wavemodesatfrequency.html) Create the wave-mode `G` EVP at fixed frequency.
  + [`waveModesAtWavenumber`](/internal-modes/classes-v2/core/imeigenvalueproblem/wavemodesatwavenumber.html) Create the wave-mode `G` EVP at fixed horizontal wavenumber.
+ Build custom EVPs
  + [`IMEigenvalueProblem`](/internal-modes/classes-v2/core/imeigenvalueproblem/imeigenvalueproblem.html) Create a physical-coordinate EVP descriptor.
+ Assemble EVPs
  + [`assemble`](/internal-modes/classes-v2/core/imeigenvalueproblem/assemble.html) Build generalized-EVP matrices on a solver's native basis.
  + [`boundaryConditions`](/internal-modes/classes-v2/core/imeigenvalueproblem/boundaryconditions.html) Placed boundary-condition array.
  + [`contextForSolver`](/internal-modes/classes-v2/core/imeigenvalueproblem/contextforsolver.html) Return the coefficient context for this EVP and solver.
  + [`leftOperator`](/internal-modes/classes-v2/core/imeigenvalueproblem/leftoperator.html) Left differential operator.
  + [`rightOperator`](/internal-modes/classes-v2/core/imeigenvalueproblem/rightoperator.html) Right differential operator.
+ Inspect EVP metadata
  + [`defaultNormalization`](/internal-modes/classes-v2/core/imeigenvalueproblem/defaultnormalization.html) Natural default normalization for this EVP.
  + [`f0`](/internal-modes/classes-v2/core/imeigenvalueproblem/f0.html) Coriolis parameter.
  + [`formulation`](/internal-modes/classes-v2/core/imeigenvalueproblem/formulation.html) Solved vertical-structure formulation.
  + [`g`](/internal-modes/classes-v2/core/imeigenvalueproblem/g.html) Gravitational acceleration.
  + [`hFromEigenvalue`](/internal-modes/classes-v2/core/imeigenvalueproblem/hfromeigenvalue.html) Equivalent-depth conversion function.
  + [`innerWeights`](/internal-modes/classes-v2/core/imeigenvalueproblem/innerweights.html) Inner-product weights for `F` and `G`.
  + [`name`](/internal-modes/classes-v2/core/imeigenvalueproblem/name.html) Short EVP name.
  + [`normalizations`](/internal-modes/classes-v2/core/imeigenvalueproblem/normalizations.html) Named modal normalization rules.
  + [`parameters`](/internal-modes/classes-v2/core/imeigenvalueproblem/parameters.html) Stored factory-specific physical inputs.
+ Select retained modes
  + [`indexValidationMode`](/internal-modes/classes-v2/core/imeigenvalueproblem/indexvalidationmode.html) Index validation behavior.
  + [`nNullModes`](/internal-modes/classes-v2/core/imeigenvalueproblem/nnullmodes.html) Number of true null modes.
  + [`partialDepthPEIndexPolicy`](/internal-modes/classes-v2/core/imeigenvalueproblem/partialdepthpeindexpolicy.html) Return the partial-depth potential-energy index policy.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Select retained modes
  + [`classifyEigenvalues`](/internal-modes/classes-v2/core/imeigenvalueproblem/classifyeigenvalues.html) Classify eigenvalues using this EVP's index metadata.
  + [`selectModes`](/internal-modes/classes-v2/core/imeigenvalueproblem/selectmodes.html) Select and label retained eigenmodes.


---