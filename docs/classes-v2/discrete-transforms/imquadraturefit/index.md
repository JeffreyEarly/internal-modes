---
layout: default
title: IMQuadratureFit
has_children: false
has_toc: false
mathjax: true
parent: Discrete transforms
grand_parent: Class documentation V2
nav_order: 2
---

#  IMQuadratureFit

Store a fixed-point quadrature fit and geometric comparison.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMQuadratureFit</code></pre></div></div>

## Overview

`IMQuadratureFit` compares two scalar Galerkin transforms on the same
physical points: one using fitted increments and one using geometric
control-volume increments. If the least-squares objective is
$$\|A\Delta z-b\|_2$$, the stored residuals are
$$A\Delta z_\mathrm{fit}-b$$ and
$$A\Delta z_\mathrm{geometric}-b$$.

Construct fits with `IMBasisSet.fitQuadrature`.

```matlab
fit = basisSet.fitQuadrature(z=z,nModes=8);
transform = fit.fittedTransform;
```




## Topics
+ Create quadrature fits
  + [`IMQuadratureFit`](/internal-modes/classes-v2/discrete-transforms/imquadraturefit/imquadraturefit.html) Create a fixed-point quadrature fit result.
+ Inspect quadrature fits
  + [`depthConstraint`](/internal-modes/classes-v2/discrete-transforms/imquadraturefit/depthconstraint.html) Whether the fit imposed exact full-depth coverage.
  + [`depthTarget`](/internal-modes/classes-v2/discrete-transforms/imquadraturefit/depthtarget.html) Requested full physical depth.
  + [`exitFlag`](/internal-modes/classes-v2/discrete-transforms/imquadraturefit/exitflag.html) Exit flag returned by `lsqlin`.
  + [`fittedDepthError`](/internal-modes/classes-v2/discrete-transforms/imquadraturefit/fitteddeptherror.html) Difference between fitted increment sum and `depthTarget`.
  + [`fittedIncrements`](/internal-modes/classes-v2/discrete-transforms/imquadraturefit/fittedincrements.html) Fitted quadrature increments.
  + [`fittedObjectiveResidual`](/internal-modes/classes-v2/discrete-transforms/imquadraturefit/fittedobjectiveresidual.html) Objective residual at the fitted increments.
  + [`fittedResidualNorm`](/internal-modes/classes-v2/discrete-transforms/imquadraturefit/fittedresidualnorm.html) Two-norm of `fittedObjectiveResidual`.
  + [`fittedTransform`](/internal-modes/classes-v2/discrete-transforms/imquadraturefit/fittedtransform.html) Transform using the fitted quadrature increments.
  + [`geometricDepthError`](/internal-modes/classes-v2/discrete-transforms/imquadraturefit/geometricdeptherror.html) Difference between geometric increment sum and `depthTarget`.
  + [`geometricIncrements`](/internal-modes/classes-v2/discrete-transforms/imquadraturefit/geometricincrements.html) Geometric control-volume increments.
  + [`geometricObjectiveResidual`](/internal-modes/classes-v2/discrete-transforms/imquadraturefit/geometricobjectiveresidual.html) Objective residual at the geometric increments.
  + [`geometricResidualNorm`](/internal-modes/classes-v2/discrete-transforms/imquadraturefit/geometricresidualnorm.html) Two-norm of `geometricObjectiveResidual`.
  + [`geometricTransform`](/internal-modes/classes-v2/discrete-transforms/imquadraturefit/geometrictransform.html) Transform using geometric control-volume increments.
  + [`nonnegativeConstraint`](/internal-modes/classes-v2/discrete-transforms/imquadraturefit/nonnegativeconstraint.html) Whether the fit imposed nonnegative increments.
  + [`objectiveMatrix`](/internal-modes/classes-v2/discrete-transforms/imquadraturefit/objectivematrix.html) Least-squares matrix $$A$$.
  + [`objectiveName`](/internal-modes/classes-v2/discrete-transforms/imquadraturefit/objectivename.html) Name of the least-squares objective.
  + [`objectiveTarget`](/internal-modes/classes-v2/discrete-transforms/imquadraturefit/objectivetarget.html) Least-squares target $$b$$.
  + [`solverOutput`](/internal-modes/classes-v2/discrete-transforms/imquadraturefit/solveroutput.html) Solver diagnostics returned by `lsqlin`.


---