---
layout: default
title: IMInternalModesQuadratureWeightFit
has_children: false
has_toc: false
mathjax: true
parent: Discrete transforms
grand_parent: Class documentation V2
nav_order: 6
---

#  IMInternalModesQuadratureWeightFit

Store diagnostics for a shared F/G quadrature-weight fit.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMInternalModesQuadratureWeightFit</code></pre></div></div>

## Overview

The fitted and geometric transforms use identical points, aligned
modal families, and requested variables.  The least-squares rows are
the unweighted stack of each variable's compressed normalized-Gram
Frobenius system. `objectiveRowVariables` records which physical
channel generated every row.
The built-in objective is

$$\sum_{V\in\mathcal V}\left\lVert
S_V(\Gamma_V(w)-\Gamma_{0,V})S_V
\right\rVert_\mathrm{F}^{2}.$$

`variableResidualNorm` and `variableGeometricResidualNorm` recover the
individual channel contributions to this stacked objective.




## Topics
+ Other
  + [`availableVariables`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/availablevariables.html) Direct channels included in the stacked fit.
  + [`depthConstraint`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/depthconstraint.html) Whether exact full-depth weight sum was imposed.
  + [`depthError`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/deptherror.html) Fitted weight-sum error relative to depth.
  + [`depthTarget`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/depthtarget.html) Full physical depth targeted by the constraint.
  + [`exitFlag`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/exitflag.html) `lsqlin` exit flag.
  + [`geometricDepthError`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/geometricdeptherror.html) Geometric weight-sum error relative to depth.
  + [`geometricResidual`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/geometricresidual.html) Geometric objective residual vector.
  + [`geometricResidualNorm`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/geometricresidualnorm.html) Norm of the geometric stacked residual.
  + [`geometricTransform`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/geometrictransform.html) Specialized transform using geometric control-volume weights.
  + [`geometricWeights`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/geometricweights.html) Geometric control-volume weights.
  + [`nonnegativeConstraint`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/nonnegativeconstraint.html) Whether nonnegative weights were imposed.
  + [`objectiveMatrix`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/objectivematrix.html) Least-squares matrix with one column per sample.
  + [`objectiveModePairs`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/objectivemodepairs.html) Aligned family-column pair for every objective row; [0 0] means custom.
  + [`objectiveName`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/objectivename.html) Name of the least-squares objective.
  + [`objectiveRowVariables`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/objectiverowvariables.html) F/G/custom provenance for every objective row.
  + [`objectiveTarget`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/objectivetarget.html) Least-squares target vector.
  + [`residual`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/residual.html) Fitted objective residual vector.
  + [`residualNorm`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/residualnorm.html) Norm of the fitted stacked residual.
  + [`solverOutput`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/solveroutput.html) `lsqlin` optimizer diagnostics.
  + [`transform`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/transform.html) Specialized transform using fitted shared weights.
  + [`variableGeometricResidualNorm`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/variablegeometricresidualnorm.html) Return the geometric normalized-Gram residual for one variable.
  + [`variableResidualNorm`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/variableresidualnorm.html) Return the fitted normalized-Gram residual for one variable.
  + [`weights`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/weights.html) Fitted shared quadrature weights.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + Construct fit results
    + [`IMInternalModesQuadratureWeightFit`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesquadratureweightfit/iminternalmodesquadratureweightfit.html) Create diagnostics for one shared quadrature-weight fit.


---