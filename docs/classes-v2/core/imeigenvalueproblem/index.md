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

Describe an internal-mode eigenvalue problem in physical coordinates.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMEigenvalueProblem</code></pre></div></div>

## Overview

`IMEigenvalueProblem` stores the physical operators, boundary laws,
formulation, inner-product weights, named normalization rules,
default normalization, equivalent-depth interpretation, and
index-selection metadata needed by coordinate-aware solvers.

```matlab
evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=1e-4);
```




## Topics
+ Create EVPs
  + [`IMEigenvalueProblem`](/internal-modes/classes-v2/core/imeigenvalueproblem/imeigenvalueproblem.html) Create a physical-coordinate EVP descriptor.
  + [`hydrostaticFModes`](/internal-modes/classes-v2/core/imeigenvalueproblem/hydrostaticfmodes.html) Create the geostrophic hydrostatic `F`-mode EVP.
  + [`hydrostaticGModes`](/internal-modes/classes-v2/core/imeigenvalueproblem/hydrostaticgmodes.html) Create the hydrostatic `G`-mode EVP.
  + [`partialDepthPEIndexPolicy`](/internal-modes/classes-v2/core/imeigenvalueproblem/partialdepthpeindexpolicy.html) Return the manuscript partial-depth PE index policy.
  + [`waveModesAtFrequency`](/internal-modes/classes-v2/core/imeigenvalueproblem/wavemodesatfrequency.html) Create the wave-mode `G` EVP at fixed frequency.
  + [`waveModesAtWavenumber`](/internal-modes/classes-v2/core/imeigenvalueproblem/wavemodesatwavenumber.html) Create the wave-mode `G` EVP at fixed horizontal wavenumber.
+ Assemble EVPs
  + [`assemble`](/internal-modes/classes-v2/core/imeigenvalueproblem/assemble.html) the EVP on a solver's native basis.
  + [`boundaryConditions`](/internal-modes/classes-v2/core/imeigenvalueproblem/boundaryconditions.html) Placed boundary-condition array.
  + [`contextForSolver`](/internal-modes/classes-v2/core/imeigenvalueproblem/contextforsolver.html) Return the coefficient context for this EVP and solver.
  + [`leftOperator`](/internal-modes/classes-v2/core/imeigenvalueproblem/leftoperator.html) Left physical operator.
  + [`rightOperator`](/internal-modes/classes-v2/core/imeigenvalueproblem/rightoperator.html) Right physical operator.
+ Inspect EVP metadata
  + [`defaultNormalization`](/internal-modes/classes-v2/core/imeigenvalueproblem/defaultnormalization.html) Natural default normalization for this EVP.
  + [`f0`](/internal-modes/classes-v2/core/imeigenvalueproblem/f0.html) Coriolis parameter owned by this EVP.
  + [`formulation`](/internal-modes/classes-v2/core/imeigenvalueproblem/formulation.html) Solved vertical-structure formulation.
  + [`g`](/internal-modes/classes-v2/core/imeigenvalueproblem/g.html) Gravitational acceleration owned by this EVP.
  + [`hFromEigenvalue`](/internal-modes/classes-v2/core/imeigenvalueproblem/hfromeigenvalue.html) Equivalent-depth conversion function.
  + [`indexValidationMode`](/internal-modes/classes-v2/core/imeigenvalueproblem/indexvalidationmode.html) Index validation behavior.
  + [`innerWeights`](/internal-modes/classes-v2/core/imeigenvalueproblem/innerweights.html) Inner-product weights for `F` and `G`.
  + [`nNullModes`](/internal-modes/classes-v2/core/imeigenvalueproblem/nnullmodes.html) Number of true null modes.
  + [`name`](/internal-modes/classes-v2/core/imeigenvalueproblem/name.html) Short EVP name.
  + [`normalizations`](/internal-modes/classes-v2/core/imeigenvalueproblem/normalizations.html) Named modal normalization rules.
  + [`parameters`](/internal-modes/classes-v2/core/imeigenvalueproblem/parameters.html) Stored factory-specific physical inputs.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Assemble EVPs
  + [`classifyEigenvalues`](/internal-modes/classes-v2/core/imeigenvalueproblem/classifyeigenvalues.html) Classify eigenvalues using this EVP's index metadata.
  + [`selectModes`](/internal-modes/classes-v2/core/imeigenvalueproblem/selectmodes.html) Select and label retained eigenmodes.


---