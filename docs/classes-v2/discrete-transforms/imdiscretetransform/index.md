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

Store one scalar discrete Galerkin transform.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMDiscreteTransform</code></pre></div></div>

## Overview

`IMDiscreteTransform` stores a sampled basis $$\Phi$$, sampled metric
$$W$$, continuous diagonal target $$\Gamma_0$$, and the Galerkin
forward matrix

$$
A_{\mathrm{gal}}=(\Phi^\mathsf{T}W\Phi)^{-1}\Phi^\mathsf{T}W.
$$

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
+ Inspect discrete transforms
  + [`basisConditionNumber`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/basisconditionnumber.html) Two-norm condition number of the sampled basis matrix.
  + [`basisMatrix`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/basismatrix.html) Sampled reconstruction basis $$\Phi$$.
  + [`forwardMatrix`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/forwardmatrix.html) Galerkin forward matrix $$A_{\mathrm{gal}}$$.
  + [`gramConditionNumber`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/gramconditionnumber.html) Two-norm condition number of the sampled Gram matrix.
  + [`gramMatrix`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/grammatrix.html) Sampled Gram matrix $$\Gamma=\Phi^\mathsf{T}W\Phi$$.
  + [`hasNegativeIncrements`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/hasnegativeincrements.html) Whether at least one supplied quadrature increment is negative.
  + [`increments`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/increments.html) Quadrature increments aligned with `z`.
  + [`metricMatrix`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/metricmatrix.html) Sample-space metric matrix $$W$$.
  + [`modeNumber`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/modenumber.html) Retained physical mode labels.
  + [`normalization`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/normalization.html) Basis-set normalization used to sample the modes.
  + [`relativeGramError`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/relativegramerror.html) Signed-norm-scaled relative Gram error.
  + [`roundTripError`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/roundtriperror.html) Coefficient round-trip error $$\|A_{\mathrm{gal}}\Phi-I\|_2$$.
  + [`targetGramIsPositiveDefinite`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/targetgramispositivedefinite.html) Whether every diagonal entry of $$\Gamma_0$$ is positive.
  + [`targetGramMatrix`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/targetgrammatrix.html) Continuous diagonal Gram target $$\Gamma_0$$.
  + [`z`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/z.html) Physical sample points.
+ Apply discrete transforms
  + [`project`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/project.html) sampled profiles onto retained modal coefficients.
  + [`reconstruct`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransform/reconstruct.html) sampled profiles from retained coefficients.


---