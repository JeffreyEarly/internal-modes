---
layout: default
title: IMQuadratureWeightFit
has_children: false
has_toc: false
mathjax: true
parent: Discrete transforms
grand_parent: Class documentation V2
nav_order: 2
---

#  IMQuadratureWeightFit

Store diagnostics for quadrature weights fitted on fixed points.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMQuadratureWeightFit</code></pre></div></div>

## Overview

`IMQuadratureWeightFit` compares a transform using fitted algebraic
weights with a transform using geometric control-volume weights on the
same fixed points. If the least-squares objective is

$$
\min_w\left\|A_{\mathrm{LS}}w-b_{\mathrm{LS}}\right\|_2,
$$

then `residual` and `geometricResidual` are

$$
A_{\mathrm{LS}}w-b_{\mathrm{LS}},
\qquad
A_{\mathrm{LS}}w_{\mathrm{geometric}}-b_{\mathrm{LS}}.
$$

Obtain the fitted weights as the primary output and this diagnostic
object as the optional second output of
`IMBasisSet.quadratureWeightsForPoints`.

```matlab
[weights,weightFit] = basisSet.quadratureWeightsForPoints(z=z,nModes=8);
transform = weightFit.transform;
```




## Topics
+ Inspect fitted weights
  + [`transform`](/internal-modes/classes-v2/discrete-transforms/imquadratureweightfit/transform.html) Discrete transform constructed with the fitted weights.
  + [`weights`](/internal-modes/classes-v2/discrete-transforms/imquadratureweightfit/weights.html) Fitted algebraic quadrature weights aligned with the fixed points.
+ Compare geometric weights
  + [`geometricTransform`](/internal-modes/classes-v2/discrete-transforms/imquadratureweightfit/geometrictransform.html) Reference transform using geometric control-volume weights.
  + [`geometricWeights`](/internal-modes/classes-v2/discrete-transforms/imquadratureweightfit/geometricweights.html) Geometric control-volume weights aligned with the fixed points.
+ Assess fit quality
  + [`objectiveName`](/internal-modes/classes-v2/discrete-transforms/imquadratureweightfit/objectivename.html) Name of the least-squares objective used to fit the weights.
  + [`residualNorm`](/internal-modes/classes-v2/discrete-transforms/imquadratureweightfit/residualnorm.html) Two-norm of the fitted objective residual.
  + [`geometricResidualNorm`](/internal-modes/classes-v2/discrete-transforms/imquadratureweightfit/geometricresidualnorm.html) Two-norm of the geometric-weight objective residual.
+ Inspect constraints
  + [`nonnegativeConstraint`](/internal-modes/classes-v2/discrete-transforms/imquadratureweightfit/nonnegativeconstraint.html) Whether the optimization required nonnegative weights.
  + [`depthConstraint`](/internal-modes/classes-v2/discrete-transforms/imquadratureweightfit/depthconstraint.html) Whether the weights were constrained to cover the full depth.
  + [`depthTarget`](/internal-modes/classes-v2/discrete-transforms/imquadratureweightfit/depthtarget.html) Full physical depth targeted by the weight-sum constraint.
  + [`depthError`](/internal-modes/classes-v2/discrete-transforms/imquadratureweightfit/deptherror.html) Difference between the fitted weight sum and `depthTarget`.
  + [`geometricDepthError`](/internal-modes/classes-v2/discrete-transforms/imquadratureweightfit/geometricdeptherror.html) Difference between the geometric weight sum and `depthTarget`.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + Construct fit results
    + [`IMQuadratureWeightFit`](/internal-modes/classes-v2/discrete-transforms/imquadratureweightfit/imquadratureweightfit.html) Create diagnostics for a quadrature-weight fit.
  + Inspect optimizer output
    + [`exitFlag`](/internal-modes/classes-v2/discrete-transforms/imquadratureweightfit/exitflag.html) Exit flag returned by `lsqlin`.
    + [`solverOutput`](/internal-modes/classes-v2/discrete-transforms/imquadratureweightfit/solveroutput.html) Solver diagnostics returned by `lsqlin`.
  + Inspect least-squares system
    + [`objectiveMatrix`](/internal-modes/classes-v2/discrete-transforms/imquadratureweightfit/objectivematrix.html) Least-squares matrix $$A_{\mathrm{LS}}$$.
    + [`objectiveTarget`](/internal-modes/classes-v2/discrete-transforms/imquadratureweightfit/objectivetarget.html) Least-squares target vector $$b_{\mathrm{LS}}$$.
    + [`residual`](/internal-modes/classes-v2/discrete-transforms/imquadratureweightfit/residual.html) Objective residual vector at the fitted weights.
    + [`geometricResidual`](/internal-modes/classes-v2/discrete-transforms/imquadratureweightfit/geometricresidual.html) Objective residual vector at the geometric weights.


---