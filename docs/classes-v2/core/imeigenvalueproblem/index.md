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

igenvalueProblem` describes the canonical Sturm-Liouville EVP

$$
ac{\partial}{\partial z}
t(p(z)\frac{\partial u}{\partial z}\right)
)u(z)=\lambda r(z)u(z).
$$

dary conditions are written in the same eigenvalue-dependent
nical form:

$$
ft[a_i u(z_i)
 p(z_i)\frac{\partial u}{\partial z}(z_i)\right]
mbda\left[c_i u(z_i)
 p(z_i)\frac{\partial u}{\partial z}(z_i)\right],
ad i\in\{\mathrm{bottom},\mathrm{surface}\}.
$$

igenvalueProblem` stores this continuous scalar problem; solvers
retize it and return an `IMBasisSet` of retained modes.

example, the Dirichlet problem

$$
ac{\partial^2 u}{\partial z^2}=\lambda u,\qquad u(-1)=u(0)=0
$$

reated with:

atlab
= IMEigenvalueProblem(zDomain=[-1 0], ...
p=@(z,~) ones(size(z)), q=@(z,~) zeros(size(z)), ...
r=@(z,~) ones(size(z)), ...
surfaceBoundary=IMBoundaryCondition.dirichlet(), ...
bottomBoundary=IMBoundaryCondition.dirichlet());
```

lver discretizes this EVP and returns the retained basis set:

atlab
er = IMSolverSpectral(nEVP=64);
sSet = solver.solveEVP(evp,nModes=4);
```




## Topics
+ Create EVPs
  + [`IMEigenvalueProblem`](/internal-modes/classes-v2/core/imeigenvalueproblem/imeigenvalueproblem.html) Create a canonical scalar EVP.
+ Define canonical coefficients
  + [`bottomBoundary`](/internal-modes/classes-v2/core/imeigenvalueproblem/bottomboundary.html) Bottom boundary condition.
  + [`p`](/internal-modes/classes-v2/core/imeigenvalueproblem/p.html) Coefficient multiplying the derivative flux.
  + [`q`](/internal-modes/classes-v2/core/imeigenvalueproblem/q.html) Coefficient multiplying the solved variable on the left side.
  + [`r`](/internal-modes/classes-v2/core/imeigenvalueproblem/r.html) Metric coefficient multiplying the eigenvalue side.
  + [`surfaceBoundary`](/internal-modes/classes-v2/core/imeigenvalueproblem/surfaceboundary.html) Surface boundary condition.
  + [`zDomain`](/internal-modes/classes-v2/core/imeigenvalueproblem/zdomain.html) Physical vertical domain.
+ Summarize EVPs
  + [`summarize`](/internal-modes/classes-v2/core/imeigenvalueproblem/summarize.html) Print a readable mathematical summary of this EVP.
+ Inspect EVP configuration
  + [`name`](/internal-modes/classes-v2/core/imeigenvalueproblem/name.html) Short EVP name.
  + [`parameters`](/internal-modes/classes-v2/core/imeigenvalueproblem/parameters.html) Additional coefficient parameters.
+ Inspect inner products
  + [`endpointWeights`](/internal-modes/classes-v2/core/imeigenvalueproblem/endpointweights.html) Return endpoint metric terms implied by active conditions.
  + [`innerProduct`](/internal-modes/classes-v2/core/imeigenvalueproblem/innerproduct.html) Return the scalar inner-product recipe.
  + [`negativeEndpointWeightCount`](/internal-modes/classes-v2/core/imeigenvalueproblem/negativeendpointweightcount.html) Count negative endpoint metric weights.
+ Inspect definiteness diagnostics
  + [`definitenessDiagnostics`](/internal-modes/classes-v2/core/imeigenvalueproblem/definitenessdiagnostics.html) Assess the norm and energy signs that control negative modes.
+ Inspect mode selection
  + [`modeSelectionDiagnostics`](/internal-modes/classes-v2/core/imeigenvalueproblem/modeselectiondiagnostics.html) Summarize negative and zero mode selection.
  + [`negativeEigenvalueBounds`](/internal-modes/classes-v2/core/imeigenvalueproblem/negativeeigenvaluebounds.html) Bound negative eigenvalues using a grid-level assessment.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`assemble`](/internal-modes/classes-v2/core/imeigenvalueproblem/assemble.html) Build the canonical matrix pair on a solver grid.
  + [`contextForSolver`](/internal-modes/classes-v2/core/imeigenvalueproblem/contextforsolver.html) Return the coefficient context for this EVP and solver.
  + [`evaluateCoefficient`](/internal-modes/classes-v2/core/imeigenvalueproblem/evaluatecoefficient.html) Evaluate a scalar, vector, or coefficient function.
  + [`makeBasisSet`](/internal-modes/classes-v2/core/imeigenvalueproblem/makebasisset.html) Create the solved scalar basis set for this EVP.
  + [`selectModes`](/internal-modes/classes-v2/core/imeigenvalueproblem/selectmodes.html) Select and label retained finite-real eigenmodes.


---