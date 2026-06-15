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

Describe a canonical scalar eigenvalue problem.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMEigenvalueProblem</code></pre></div></div>

## Overview

`IMEigenvalueProblem` is the solver-independent description of the
scalar problem
$$-(p u')' + q u = \lambda r u,$$
with endpoint conditions
$$-[a_i u-b_i(pu')]=\lambda[c_i u-d_i(pu')].$$
The EVP owns the continuous problem: the physical interval, coefficient
functions, endpoint conditions, normalization rules, and diagnostic
definiteness checks. A solver owns only the numerical choices used to
discretize this problem.

```matlab
evp = IMEigenvalueProblem(zDomain=[-1 0], ...
    p=@(z,~) ones(size(z)), q=@(z,~) zeros(size(z)), ...
    r=@(z,~) ones(size(z)), ...
    surfaceBoundary=IMBoundaryCondition.dirichlet(), ...
    bottomBoundary=IMBoundaryCondition.dirichlet());
solver = IMSolverSpectral(nEVP=64);
basisSet = solver.solveEVP(evp,nModes=4);
```




## Topics
+ Create EVPs
  + [`IMEigenvalueProblem`](/internal-modes/classes-v2/core/imeigenvalueproblem/imeigenvalueproblem.html) Create a canonical scalar EVP.
+ Define canonical coefficients
  + [`bottomBoundary`](/internal-modes/classes-v2/core/imeigenvalueproblem/bottomboundary.html) Bottom endpoint condition.
  + [`p`](/internal-modes/classes-v2/core/imeigenvalueproblem/p.html) Coefficient multiplying the derivative flux.
  + [`q`](/internal-modes/classes-v2/core/imeigenvalueproblem/q.html) Coefficient multiplying the solved variable on the left side.
  + [`r`](/internal-modes/classes-v2/core/imeigenvalueproblem/r.html) Metric coefficient multiplying the eigenvalue side.
  + [`surfaceBoundary`](/internal-modes/classes-v2/core/imeigenvalueproblem/surfaceboundary.html) Surface endpoint condition.
  + [`zDomain`](/internal-modes/classes-v2/core/imeigenvalueproblem/zdomain.html) Physical vertical domain.
+ Inspect EVP configuration
  + [`defaultNormalization`](/internal-modes/classes-v2/core/imeigenvalueproblem/defaultnormalization.html) Natural default normalization.
  + [`name`](/internal-modes/classes-v2/core/imeigenvalueproblem/name.html) Short EVP name.
  + [`normalizations`](/internal-modes/classes-v2/core/imeigenvalueproblem/normalizations.html) Named normalization rules.
  + [`parameters`](/internal-modes/classes-v2/core/imeigenvalueproblem/parameters.html) Additional coefficient parameters.
+ Inspect inner products
  + [`endpointWeights`](/internal-modes/classes-v2/core/imeigenvalueproblem/endpointweights.html) Return endpoint metric terms implied by active conditions.
  + [`innerProduct`](/internal-modes/classes-v2/core/imeigenvalueproblem/innerproduct.html) Return the scalar inner-product recipe.
  + [`negativeEndpointWeightCount`](/internal-modes/classes-v2/core/imeigenvalueproblem/negativeendpointweightcount.html) Count negative endpoint metric weights.
+ Inspect definiteness diagnostics
  + [`definitenessInfo`](/internal-modes/classes-v2/core/imeigenvalueproblem/definitenessinfo.html) Check grid-level signs for the canonical coefficients.
+ Inspect mode selection
  + [`modeSelectionDiagnostics`](/internal-modes/classes-v2/core/imeigenvalueproblem/modeselectiondiagnostics.html) Summarize negative and zero mode selection.
  + [`negativeEigenvalueBounds`](/internal-modes/classes-v2/core/imeigenvalueproblem/negativeeigenvaluebounds.html) Bound negative eigenvalues using a grid-level assessment.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`assemble`](/internal-modes/classes-v2/core/imeigenvalueproblem/assemble.html) Build the canonical matrix pair on a solver grid.
  + [`contextForSolver`](/internal-modes/classes-v2/core/imeigenvalueproblem/contextforsolver.html) Return the coefficient context for this EVP and solver.
  + [`coordinateProfile`](/internal-modes/classes-v2/core/imeigenvalueproblem/coordinateprofile.html) Return fields needed by a solver coordinate map.
  + [`evaluateCoefficient`](/internal-modes/classes-v2/core/imeigenvalueproblem/evaluatecoefficient.html) Evaluate a scalar, vector, or coefficient function.
  + [`makeBasisSet`](/internal-modes/classes-v2/core/imeigenvalueproblem/makebasisset.html) Create the solved scalar basis set for this EVP.
  + [`selectModes`](/internal-modes/classes-v2/core/imeigenvalueproblem/selectmodes.html) Select and label retained finite-real eigenmodes.


---