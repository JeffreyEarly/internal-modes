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

A quadrature-weight fit holds the physical sample points $$z_k$$, the
retained modes, and their normalization fixed. Only one quadrature
weight $$w_k$$ per point is optimized. For the sampled mode matrix

$$
\Phi_{ki}=u_i(z_k),
$$

the weights produce the discrete Gram matrix

$$
\Gamma(w)
=\Phi^{\mathsf T}
\operatorname{diag}\!\left(r(z_k)w_k\right)\Phi
+\Gamma_{\mathrm{endpoint}}.
$$

Here $$r(z)$$ is the EVP interior inner-product weight and
$$\Gamma_{\mathrm{endpoint}}$$ contains endpoint contributions that do
not depend on the interior quadrature weights. The continuous Gram
target $$\Gamma_0$$ contains the desired inner products of the retained
normalized modes. The default objective adjusts $$w$$ so that
$$\Gamma(w)\approx\Gamma_0$$. With
$$C_i=(\Gamma_0)_{ii}$$, its least-squares system is

$$
(A_{\mathrm{LS}})_{(i,j),k}
=\rho_{ij}\frac{r(z_k)\Phi_{ki}\Phi_{kj}}{\sqrt{|C_iC_j|}},
\qquad
(b_{\mathrm{LS}})_{(i,j)}
=\rho_{ij}\frac{(\Gamma_0-\Gamma_{\mathrm{endpoint}})_{ij}}
{\sqrt{|C_iC_j|}},
$$

for $$1\leq i\leq j\leq n_m$$, where

$$
\rho_{ij}=
\begin{cases}
1, & i=j,\\
\sqrt{2}, & i<j.
\end{cases}
$$

The fitted weights solve

$$
\min_w\left\|A_{\mathrm{LS}}w-b_{\mathrm{LS}}\right\|_2.
$$

Equivalently, define

$$
S=\operatorname{diag}\!\left(
\left|\operatorname{diag}\Gamma_0\right|^{-1/2}
\right),
\qquad
E(w)=S\left(\Gamma(w)-\Gamma_0\right)S.
$$

The default `"normalizedGramFrobenius"` objective minimizes
$$\|E(w)\|_{\mathrm F}$$. This is an aggregate error over the full
retained Gram matrix: diagonal entries measure errors in individual
modal norms, and off-diagonal entries measure lost orthogonality. The
least-squares system stores only the upper triangle. Its $$\sqrt{2}$$
off-diagonal factor preserves the exact identity

$$
\|E(w)\|_{\mathrm F}^2
=\sum_i E_{ii}(w)^2+2\sum_{i<j}E_{ij}(w)^2.
$$

By default the fit also requires

$$
w_k\geq0,
\qquad
\sum_k w_k=D,
$$

where $$D=z_s-z_b$$ is the full physical depth.

Before optimization, the fixed points define a natural geometric
control-volume rule. For $$n$$ points, let

$$
e_1=z_b,
\qquad
e_{k+1}=\frac{z_k+z_{k+1}}{2}\quad(k=1,\ldots,n-1),
\qquad
e_{n+1}=z_s.
$$

The unoptimized weights are

$$
w_k^{\mathrm{geometric}}=e_{k+1}-e_k.
$$

These are the literal midpoint/control-volume approximation to the
sampled integral: they use only point geometry and do not use modal
orthogonality. `quadratureWeightsForPoints` supplies them as the
optimizer's initial guess. They are positive, cover the full domain,
and define `geometricTransform`, the reference obtained by skipping
optimization. The fitted algebraic weights are instead adjusted to
reproduce the retained modal inner products and need not retain a
geometric-width interpretation when constraints are relaxed.

`transform` and `weights` are the optimized production result;
`geometricTransform` and `geometricWeights` are the unoptimized
comparison. For the built-in objective, `residualNorm` and
`geometricResidualNorm` are the Frobenius norms of $$E$$ for those two
rules. For a custom objective they instead mean the generic residual
norm $$\|Aw-b\|_2$$. Their transforms'
`relativeGramOperatorError` values report $$\|E\|_2$$, the largest
Gram distortion over any normalized combination of retained modes.
`roundTripError` is different again: it measures algebraic recovery
of retained coefficients and can be tiny even when the sampled
quadrature does not accurately reproduce the continuous Gram matrix.
Constraint properties record whether nonnegativity and full-depth
coverage were imposed, while `depthError` and
`transform.hasNegativeWeights` report the corresponding fitted result.
A custom objective changes $$A_{\mathrm{LS}}$$ and
$$b_{\mathrm{LS}}$$, but the fitted and geometric rules still use the
same points, modes, and normalization. The quadrature-weight
regression sweep supports retaining the unregularized Frobenius
objective as the default; geometric weights remain the optimizer's
initial guess and the unoptimized comparison baseline.

Obtain the fitted weights as the primary output and this diagnostic
object as the optional second output of
`IMBasisSet.quadratureWeightsForPoints`. The returned `weights`,
`weightFit.weights`, and `weightFit.transform.weights` are the same
vector.

```matlab
[weights,weightFit] = basisSet.quadratureWeightsForPoints(z=z,nModes=8);
weightFit.residualNorm
weightFit.geometricResidualNorm
weightFit.transform.relativeGramOperatorError
weightFit.geometricTransform.relativeGramOperatorError
weightFit.transform.roundTripError
coefficients = weightFit.transform.transformForward(values);
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
  + [`residualNorm`](/internal-modes/classes-v2/discrete-transforms/imquadratureweightfit/residualnorm.html) Norm of the fitted objective residual.
  + [`geometricResidualNorm`](/internal-modes/classes-v2/discrete-transforms/imquadratureweightfit/geometricresidualnorm.html) Norm of the geometric-weight objective residual.
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