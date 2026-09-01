---
layout: default
title: IMDiscreteTransformAssessment
has_children: false
has_toc: false
mathjax: true
parent: Discrete transforms
grand_parent: Class documentation V2
nav_order: 2
---

#  IMDiscreteTransformAssessment

Store fixed-rule retained-band diagnostics for a scalar transform.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMDiscreteTransformAssessment</code></pre></div></div>

## Overview

`IMDiscreteTransformAssessment` records one final point-and-weight rule
and the leading modal prefixes assessed on that rule. For certified
construction, `certificationSearch` records the independently refitted
counts considered before this exact final band was selected, and
`weightFitModeCount` equals `retainedModeCount`. Lower-level fixed-rule
calls can still fit a larger candidate once and return a shorter
prefix; their differing counts make that behavior explicit.

For prefix $$N$$, let $$\Gamma_N$$ be its sampled Gram matrix,
$$\Gamma_{0,N}$$ its diagonal continuous target, and

$$
S_N=\operatorname{diag}\!\left(
\left|\operatorname{diag}\Gamma_{0,N}\right|^{-1/2}
\right),\qquad
E_N=S_N(\Gamma_N-\Gamma_{0,N})S_N.
$$

The default Gram policy uses
$$\lVert E_N\rVert_2$$, stored in `gramError`. By contrast, the
built-in quadrature objective reports $$\lVert E_N\rVert_\mathrm F$$
as `IMQuadratureWeightFit.residualNorm`. `roundTripError` is
$$\lVert A_\mathrm fA_\mathrm i-I\rVert_2$$; it measures algebraic
consistency and can be near roundoff even when the Gram rule is poor.
The two condition-number columns describe algebraic sensitivity and
are diagnostics rather than acceptance policies.

When enabled, rejected-mode leakage is

$$
\ell_N=\max_{N<j\leq N_\mathrm{check}}
\frac{\lVert\Pi_N^\mathrm{discrete}u_j\rVert_\mu}
{\lVert u_j\rVert_\mu},
$$

and scalar quadratic aliasing is

$$
q_N=\max_{1\leq i\leq j\leq N}
\frac{\left\lVert
\Pi_N^\mathrm{discrete}(u_i u_j)-
\Pi_N^\mathrm{continuous}(u_i u_j)
\right\rVert_\mu}{\lVert u_i u_j\rVert_\mu}.
$$

The continuous quadratic reference is integrated on the source
solver's inner-product grid, independently of the fitted points and
weights. Leakage and quadratic aliasing require a positive-definite
target because $$\lVert\cdot\rVert_\mu$$ must be a true norm. Signed
targets retain Gram diagnostics but cannot enable those policies.

`prefixDiagnostics` has one row per candidate prefix. Its columns are:

- `modeCount`: prefix column count $$N$$.
- `lastModeNumber`: physical label of column $$N$$.
- `gramError`: $$\lVert E_N\rVert_2$$.
- `roundTripError`: $$\lVert A_\mathrm fA_\mathrm i-I\rVert_2$$.
- `inverseMatrixConditionNumber`: $$\kappa_2(A_\mathrm i)$$.
- `gramConditionNumber`: $$\kappa_2(\Gamma_N)$$.
- `sampledGramRank`: numerical rank of $$\Gamma_N$$.
- `leakageError`: $$\ell_N$$, or `NaN` when disabled.
- `leakageLimitingModeNumber`: rejected physical mode label attaining $$\ell_N$$.
- `quadraticAliasingError`: $$q_N$$, or `NaN` when disabled.
- `quadraticLimitingModeNumberI`, `quadraticLimitingModeNumberJ`: physical labels of the product attaining $$q_N$$.
- `gramAccepted`, `leakageAccepted`, `quadraticAccepted`: cumulative per-policy prefix decisions.
- `combinedAccepted`: intersection of every enabled cumulative policy.

Acceptance is consecutive: after a policy first rejects a prefix,
larger prefixes remain rejected even if a later raw metric happens to
fall below tolerance. Each policy struct stores `enabled`, `tolerance`,
the per-prefix `error` and `accepted` arrays, `maximumAcceptedModeCount`,
`limitingValue`, and a readable `reason`. Leakage additionally stores
`nCheckModes` and its limiting mode labels; quadratic aliasing stores
its limiting product labels.

```matlab
[transform,assessment] = basisSet.discreteTransform(nPoints=32);
assessment.prefixDiagnostics
assessment.gramPolicy
assessment.leakagePolicy
assessment.quadraticAliasingPolicy
```




## Topics
+ Create transform assessments
  + [`IMDiscreteTransformAssessment`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransformassessment/imdiscretetransformassessment.html) Create a scalar discrete-transform assessment.
+ Inspect transform assessments
  + [`actualPointCount`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransformassessment/actualpointcount.html) Actual physical sample count.
  + [`candidateModeCount`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransformassessment/candidatemodecount.html) Number of modes in the full candidate band.
  + [`candidateModeNumber`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransformassessment/candidatemodenumber.html) Physical labels for the candidate modes.
  + [`candidateTransform`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransformassessment/candidatetransform.html) Full candidate-band transform built with the fixed rule.
  + [`prefixDiagnostics`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransformassessment/prefixdiagnostics.html) Per-prefix transform and policy diagnostics.
  + [`requestedPointCount`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransformassessment/requestedpointcount.html) Requested physical sample count.
  + [`retainedModeCount`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransformassessment/retainedmodecount.html) Number of modes retained by the combined policies.
  + [`retainedModeNumber`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransformassessment/retainedmodenumber.html) Physical labels for the retained modes.
  + [`transform`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransformassessment/transform.html) Production transform using the retained modal prefix.
  + [`weightFit`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransformassessment/weightfit.html) Quadrature-weight fit used to build the fixed rule.
+ Inspect retained-band policies
  + [`gramPolicy`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransformassessment/grampolicy.html) Normalized-Gram retained-band policy result.
  + [`leakagePolicy`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransformassessment/leakagepolicy.html) Rejected-mode leakage policy result.
  + [`limitingPolicy`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransformassessment/limitingpolicy.html) Policy that imposed the final retained count.
  + [`quadraticAliasingPolicy`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransformassessment/quadraticaliasingpolicy.html) Scalar quadratic-aliasing policy result.
  + [`retentionReason`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransformassessment/retentionreason.html) Readable explanation of the retained-band decision.
+ Other
  + [`certificationSearch`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransformassessment/certificationsearch.html) Independently refitted count-selection attempts.
  + [`gridDesign`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransformassessment/griddesign.html) Provenance of the physical sample grid.
  + [`weightFitModeCount`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransformassessment/weightfitmodecount.html) Number of modes used when fitting the stored quadrature weights.
  + [`withCertificationMetadata`](/internal-modes/classes-v2/discrete-transforms/imdiscretetransformassessment/withcertificationmetadata.html) Return this assessment with grid and count-search provenance.


---