---
layout: default
title: IMInternalModesDiscreteTransform
has_children: false
has_toc: false
mathjax: true
parent: Discrete transforms
grand_parent: Class documentation V2
nav_order: 4
---

#  IMInternalModesDiscreteTransform

Store aligned discrete transforms for internal-mode F/G variables.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMInternalModesDiscreteTransform</code></pre></div></div>

## Overview

A single instance stores one point-and-weight rule and the complete
aligned internal-mode family.  The sampled `F` and `G` columns share
physical mode labels, equivalent depths, and normalization, while
each directly representable variable owns its metric, Gram matrix,
active-column projector, and Galerkin forward matrix.

Use variable-qualified accessors because the two channels generally
have different continuous inner products:

Let $$n_z$$ be the sample count and $$n_m$$ the aligned family count.
`inverseMatrix(variable=V)` is $$n_z\times n_m$$, while
`forwardMatrix(variable=V)` is $$n_m\times n_z$$ and has zero rows
for inactive variable columns. For the active projector $$Q_V$$,

$$A_\mathrm{f}^{V}A_\mathrm{i}^{V}=Q_V.$$

`modeProjectionFunctional` returns
$$(A_\mathrm{i}^{V})^\mathsf{T}W_VX$$ before the active Gram solve;
`transformForward` returns coefficients after that solve. Both use
the signed Pontryagin pairing. `targetMajorantGramMatrix` returns the
induced positive Hilbert majorant used for error magnitudes; it does
not replace the signed projection metric.
When `variable` is omitted, accessors use `primaryVariable`, which is
the solved formulation when that channel was requested and otherwise
the first directly representable requested channel.

```matlab
[transform,assessment] = basisSet.discreteTransform(nPoints=24,variables=["F","G"]);
aG = transform.transformForward(G,variable="G");
F = transform.transformBack(aG,variable="F");
```




## Topics
+ Inspect samples and modes
  + [`primaryVariable`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/primaryvariable.html) Default F/G channel used when `variable` is omitted.
+ Assess transform quality
  + [`targetMajorantGramMatrix`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/targetmajorantgrammatrix.html) Return the continuous positive Hilbert-majorant Gram matrix.
+ Other
  + [`N2Values`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/n2values.html) Buoyancy frequency squared sampled at `z`.
  + [`activeModeMask`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/activemodemask.html) Return the full-family active-column projector mask.
  + [`availableVariables`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/availablevariables.html) Directly representable forward channels in canonical order.
  + [`depth`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/depth.html) Full physical depth.
  + [`endpointLocations`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/endpointlocations.html) Canonical endpoint identities, `["surface";"bottom"]`.
  + [`endpointValues`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/endpointvalues.html) Return normalized endpoint traces, surface then bottom.
  + [`forwardMatrix`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/forwardmatrix.html) Return the variable-qualified Galerkin forward matrix.
  + [`forwardTransformReason`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/forwardtransformreason.html) Return the reason a direct projection is unavailable.
  + [`g`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/g.html) Gravitational acceleration.
  + [`gramConditionNumber`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/gramconditionnumber.html) Return the active sampled-Gram condition number.
  + [`gramMatrix`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/grammatrix.html) Return the sampled signed full-family Gram matrix.
  + [`h`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/h.html) Equivalent depths aligned with `modeNumber`.
  + [`hasForwardTransform`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/hasforwardtransform.html) Return whether a direct sampled projection exists.
  + [`hasNegativeWeights`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/hasnegativeweights.html) Whether any shared quadrature weight is negative.
  + [`inverseMatrix`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/inversematrix.html) Return the sampled modal synthesis matrix for F or G.
  + [`inverseMatrixConditionNumber`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/inversematrixconditionnumber.html) Return the active sampled-basis condition number.
  + [`metricMatrix`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/metricmatrix.html) Return the variable-qualified sampled signed metric.
  + [`modeFamily`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/modefamily.html) Internal-mode family name.
  + [`modeNumber`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/modenumber.html) Physical labels for the aligned modal family.
  + [`modeProjectionFunctional`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/modeprojectionfunctional.html) Apply `(A_i^V)' W_V` without solving the modal Gram system.
  + [`normalization`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/normalization.html) captured when the transform was built.
  + [`problemMetadata`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/problemmetadata.html) Immutable basis metadata snapshot.
  + [`relativeGramOperatorError`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/relativegramoperatorerror.html) Return the normalized Gram operator error for one channel.
  + [`roundTripError`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/roundtriperror.html) Return the active-projector round-trip error for one channel.
  + [`targetGramIsPositiveDefinite`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/targetgramispositivedefinite.html) Return whether the active target channel defines a norm.
  + [`targetGramMatrix`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/targetgrammatrix.html) Return the continuous signed full-family target Gram matrix.
  + [`transformBack`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/transformback.html) Synthesize F or G values from full-family coefficients.
  + [`transformForward`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/transformforward.html) Project sampled values into aligned family coefficients.
  + [`weights`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/weights.html) Shared quadrature weights.
  + [`z`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/z.html) Shared physical sample locations.
  + [`zDomain`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/zdomain.html) Physical vertical domain `[zBottom zSurface]`.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + Construct transforms
    + [`IMInternalModesDiscreteTransform`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/iminternalmodesdiscretetransform.html) Create an aligned internal-mode transform from prepared data.
+ Other
  + [`channelDiagnostics`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransform/channeldiagnostics.html) Return scalar quality diagnostics for one built channel.


---