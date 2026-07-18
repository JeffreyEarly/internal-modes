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

For sampled values $$\mathbf{x}$$ and retained modal coefficients
$$\mathbf{a}$$, the transform matrices satisfy

$$
\mathbf{a}=A_{\mathrm f}\mathbf{x},
\qquad
\widehat{\mathbf{x}}=A_{\mathrm i}\mathbf{a},
$$

where `forwardMatrix` is $$A_{\mathrm f}$$ and `inverseMatrix` is
$$A_{\mathrm i}=\Phi$$, the sampled modal basis. The Galerkin forward
matrix is

$$
A_{\mathrm f}=(\Phi^\mathsf{T}W\Phi)^{-1}\Phi^\mathsf{T}W.
$$

The matrices may be rectangular. They obey

$$
A_{\mathrm f}A_{\mathrm i}=I,
\qquad
A_{\mathrm i}A_{\mathrm f}=P_W,
$$

where $$P_W$$ is the $$W$$-Galerkin projector onto the retained
sampled modal subspace. When $$W$$ is positive definite, this is the
$$W$$-orthogonal projector.

Construct transforms from solved scalar modes with
`IMBasisSet.discreteTransform`.

```matlab
transform = basisSet.discreteTransform(z=z,nModes=8);
coefficients = transform.project(values);
valuesFit = transform.reconstruct(coefficients);
```




## Topics
+ Create discrete transforms
  + [`IMDiscreteTransform`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/imdiscretetransform.html) Create a scalar discrete transform from canonical matrices.
+ Use transform matrices
  + [`forwardMatrix`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/forwardmatrix.html) Forward transform matrix $$A_{\mathrm f}$$.
  + [`inverseMatrix`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/inversematrix.html) Inverse transform matrix $$A_{\mathrm i}=\Phi$$.
+ Apply discrete transforms
  + [`project`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/project.html) Return retained modal coefficients for sampled profiles.
  + [`reconstruct`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/reconstruct.html) Return sampled profiles reconstructed from retained coefficients.
+ Inspect samples and modes
  + [`z`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/z.html) Physical sample points.
  + [`increments`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/increments.html) Quadrature increments aligned with `z`.
  + [`modeNumber`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/modenumber.html) Retained physical mode labels.
  + [`normalization`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/normalization.html) Basis-set normalization used to sample the modes.
+ Assess transform quality
  + [`relativeGramError`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/relativegramerror.html) Signed-norm-scaled relative Gram error.
  + [`roundTripError`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/roundtriperror.html) Coefficient round-trip error $$\|A_{\mathrm f}A_{\mathrm i}-I\|_2$$.
  + [`inverseMatrixConditionNumber`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/inversematrixconditionnumber.html) Two-norm condition number of `inverseMatrix`.
  + [`gramConditionNumber`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/gramconditionnumber.html) Two-norm condition number of the sampled Gram matrix.
  + [`targetGramIsPositiveDefinite`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/targetgramispositivedefinite.html) Whether every diagonal entry of $$\Gamma_0$$ is positive.
  + [`hasNegativeIncrements`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/hasnegativeincrements.html) Whether at least one supplied quadrature increment is negative.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + Inspect metric construction
    + [`metricMatrix`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/metricmatrix.html) Sample-space metric matrix $$W$$.
    + [`gramMatrix`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/grammatrix.html) Sampled Gram matrix $$\Gamma=A_{\mathrm i}^\mathsf{T}WA_{\mathrm i}$$.
    + [`targetGramMatrix`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/targetgrammatrix.html) Continuous diagonal Gram target $$\Gamma_0$$.


---