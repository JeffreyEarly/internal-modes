---
layout: default
title: IMDiscreteTransform
has_children: false
has_toc: false
mathjax: true
parent: Discrete transforms
grand_parent: Class documentation V2
nav_order: 1
---

#  IMDiscreteTransform

Store forward and inverse matrices for a scalar modal transform.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMDiscreteTransform</code></pre></div></div>

## Overview

Let $$n_z$$ be the number of sample points, $$n_m$$ the number of
retained modes, and $$n_p$$ the number of sampled profiles. The
sampled modal basis is

$$
A_{\mathrm i}=\Phi\in\mathbb{R}^{n_z\times n_m},
\qquad
\Phi_{ij}=u_j(z_i),
$$

where the columns contain the retained, normalized modes. For the
sample-space metric $$W$$, the sampled Gram matrix and Galerkin
forward matrix are

$$
\Gamma=\Phi^\mathsf{T}W\Phi,
\qquad
A_{\mathrm f}=\Gamma^{-1}\Phi^\mathsf{T}W
=\left(A_{\mathrm i}^\mathsf{T}WA_{\mathrm i}\right)^{-1}
A_{\mathrm i}^\mathsf{T}W.
$$

For sampled profiles $$X\in\mathbb{R}^{n_z\times n_p}$$ and modal
coefficients $$A\in\mathbb{R}^{n_m\times n_p}$$, the forward and
back transforms are

$$
A=A_{\mathrm f}X,
\qquad
\widehat{X}=A_{\mathrm i}A.
$$

The matrices are generally rectangular. For a full-rank retained
basis they obey

$$
A_{\mathrm f}A_{\mathrm i}=I,
\qquad
A_{\mathrm i}A_{\mathrm f}=P_W,
$$

where $$P_W$$ is the sampled-space Galerkin projector onto the
retained modal subspace. Thus transforming retained coefficients back
and then forward recovers those coefficients exactly, while a general
sampled profile is projected rather than necessarily reproduced. When
$$W$$ is positive definite, $$P_W$$ is the $$W$$-orthogonal projector.

Transform construction and quadrature accuracy are related but
distinct. If $$\Gamma_0$$ is the continuous target Gram matrix, define

$$
S=\operatorname{diag}\!\left(
\left|\operatorname{diag}\Gamma_0\right|^{-1/2}
\right),
\qquad
E=S(\Gamma-\Gamma_0)S.
$$

`relativeGramOperatorError` is $$\|E\|_2$$: the largest Gram
distortion over any normalized combination of retained modes. By
contrast, `roundTripError` measures only how accurately the two
transform matrices recover retained coefficients. It can be near
roundoff even when the quadrature reproduces $$\Gamma_0$$ poorly.

Construct transforms from solved scalar modes with
`IMBasisSet.discreteTransform`.

```matlab
transform = basisSet.discreteTransform(z=z,nModes=8);
coefficients = transform.transformForward(values);
valuesFit = transform.transformBack(coefficients);
```




## Topics
+ Create discrete transforms
  + [`IMDiscreteTransform`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/imdiscretetransform.html) Create a scalar discrete transform from canonical matrices.
+ Use transform matrices
  + [`forwardMatrix`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/forwardmatrix.html) Map sampled profiles to retained modal coefficients.
  + [`inverseMatrix`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/inversematrix.html) Map retained modal coefficients back to sampled profiles.
+ Apply discrete transforms
  + [`transformForward`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/transformforward.html) Transform sampled profiles forward to modal coefficients.
  + [`transformBack`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/transformback.html) Transform modal coefficients back to sampled profiles.
+ Inspect samples and modes
  + [`z`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/z.html) Physical points at which profiles and modes are sampled.
  + [`weights`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/weights.html) Quadrature weights associated with the sample points.
  + [`modeNumber`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/modenumber.html) Physical labels for the retained modal rows and columns.
  + [`normalization`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/normalization.html) Name of the normalization captured by this transform.
+ Assess transform quality
  + [`relativeGramOperatorError`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/relativegramoperatorerror.html) Measure the worst normalized Gram distortion.
  + [`roundTripError`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/roundtriperror.html) Measure recovery of retained modal coefficients.
  + [`inverseMatrixConditionNumber`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/inversematrixconditionnumber.html) Two-norm condition number of the sampled modal basis.
  + [`sampledGramRank`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/sampledgramrank.html) Numerical rank of the sampled modal Gram matrix.
  + [`gramConditionNumber`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/gramconditionnumber.html) Two-norm condition number of the sampled Gram matrix.
  + [`targetGramIsPositiveDefinite`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/targetgramispositivedefinite.html) Whether the target modal Gram matrix defines a positive norm.
  + [`hasNegativeWeights`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/hasnegativeweights.html) Whether at least one quadrature weight is negative.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + Inspect metric construction
    + [`metricMatrix`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/metricmatrix.html) Sample-space bilinear-form matrix $$W$$.
    + [`gramMatrix`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/grammatrix.html) Sampled modal Gram matrix $$\Gamma$$.
    + [`targetGramMatrix`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/targetgrammatrix.html) Continuous diagonal Gram matrix targeted by the quadrature rule.


---