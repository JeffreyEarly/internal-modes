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
The EVP owns the physical interval, coefficient functions, endpoint
conditions, equivalent-depth mapping, normalization rules, and diagnostic
definiteness checks.
A solver owns the numerical grid, coordinate mapping, and derivative
matrices used to discretize this continuous problem.




## Topics
+ Create EVPs
  + [`IMEigenvalueProblem`](/internal-modes/classes-v2/core/imeigenvalueproblem/imeigenvalueproblem.html) Create a canonical scalar EVP.
+ Assemble EVPs
  + [`assemble`](/internal-modes/classes-v2/core/imeigenvalueproblem/assemble.html) the canonical matrix pair on a solver grid.
  + [`bottomBoundary`](/internal-modes/classes-v2/core/imeigenvalueproblem/bottomboundary.html) Bottom endpoint condition.
  + [`contextForSolver`](/internal-modes/classes-v2/core/imeigenvalueproblem/contextforsolver.html) Return the coefficient context for this EVP and solver.
  + [`p`](/internal-modes/classes-v2/core/imeigenvalueproblem/p.html) Coefficient multiplying the derivative flux.
  + [`q`](/internal-modes/classes-v2/core/imeigenvalueproblem/q.html) Coefficient multiplying the solved variable on the left side.
  + [`r`](/internal-modes/classes-v2/core/imeigenvalueproblem/r.html) Metric coefficient multiplying the eigenvalue side.
  + [`surfaceBoundary`](/internal-modes/classes-v2/core/imeigenvalueproblem/surfaceboundary.html) Surface endpoint condition.
  + [`zDomain`](/internal-modes/classes-v2/core/imeigenvalueproblem/zdomain.html) Physical vertical domain.
+ Inspect diagnostics
  + [`defaultNormalization`](/internal-modes/classes-v2/core/imeigenvalueproblem/defaultnormalization.html) Natural default normalization.
  + [`definitenessInfo`](/internal-modes/classes-v2/core/imeigenvalueproblem/definitenessinfo.html) Check grid-level signs for the canonical coefficients.
  + [`endpointWeights`](/internal-modes/classes-v2/core/imeigenvalueproblem/endpointweights.html) Return endpoint metric terms implied by active conditions.
  + [`hFromEigenvalue`](/internal-modes/classes-v2/core/imeigenvalueproblem/hfromeigenvalue.html) Equivalent-depth conversion function.
  + [`innerProduct`](/internal-modes/classes-v2/core/imeigenvalueproblem/innerproduct.html) Return the scalar inner-product recipe.
  + [`metadata`](/internal-modes/classes-v2/core/imeigenvalueproblem/metadata.html) Additional EVP metadata.
  + [`metricIndex`](/internal-modes/classes-v2/core/imeigenvalueproblem/metricindex.html) Count negative endpoint directions in the metric.
  + [`name`](/internal-modes/classes-v2/core/imeigenvalueproblem/name.html) Short EVP name.
  + [`negativeEigenvalueBounds`](/internal-modes/classes-v2/core/imeigenvalueproblem/negativeeigenvaluebounds.html) Bound negative eigenvalues using grid-level certification.
  + [`normalizations`](/internal-modes/classes-v2/core/imeigenvalueproblem/normalizations.html) Named normalization rules.
+ Select modes
  + [`hasZeroMode`](/internal-modes/classes-v2/core/imeigenvalueproblem/haszeromode.html) Whether the scalar problem declares a zero mode.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Assemble EVPs
  + [`coordinateProfile`](/internal-modes/classes-v2/core/imeigenvalueproblem/coordinateprofile.html) Return fields needed by a solver coordinate map.
+ Select modes
  + [`selectModes`](/internal-modes/classes-v2/core/imeigenvalueproblem/selectmodes.html) Select and label retained finite-real eigenmodes.
+ Developer topics
  + [`evaluateCoefficient`](/internal-modes/classes-v2/core/imeigenvalueproblem/evaluatecoefficient.html) Evaluate a scalar, vector, or coefficient function.
  + [`makeBasisSet`](/internal-modes/classes-v2/core/imeigenvalueproblem/makebasisset.html) Create the solved scalar basis set for this EVP.


---