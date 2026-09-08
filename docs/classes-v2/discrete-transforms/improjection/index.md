---
layout: default
title: IMProjection
has_children: false
has_toc: false
mathjax: true
parent: Discrete transforms
grand_parent: Class documentation V2
nav_order: 3
---

#  IMProjection

Construct a fixed sampled projection and measure its numerical quality.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMProjection</code></pre></div></div>

## Overview

The signed sampled pairing defines coefficients through
$$A = (\Phi^T W\Phi)^{-1}\Phi^T W.$$
Rank-deficient pairings use a pseudoinverse and report infinite Gram
error. The separately supplied positive majorant measures coefficient
errors; it never replaces the physical pairing. This object does not
choose a grid, fit weights, accept tolerances, or select model counts.
`fromPairing` instead accepts an explicit physical source dual and its
coefficient system. It has no synthesis basis and cannot measure a
sampled-basis Gram discrepancy or round trip.

```matlab
projection = IMProjection(sampledBasis,metricMatrix,targetGramMatrix);
coefficients = projection.project(values);
```




## Topics
+ Create projections
  + [`IMProjection`](/internal-modes/classes-v2/discrete-transforms/improjection/improjection.html) Prepare a projection without making acceptance decisions.
  + [`fromPairing`](/internal-modes/classes-v2/discrete-transforms/improjection/frompairing.html) Construct a prescribed physical source dual without inventing synthesis data.
  + [`prefix`](/internal-modes/classes-v2/discrete-transforms/improjection/prefix.html) Construct the projection of a leading set of array columns.
+ Inspect projection data
  + [`activeColumnMask`](/internal-modes/classes-v2/discrete-transforms/improjection/activecolumnmask.html) Logical row identifying columns with a direct projection.
  + [`columnCount`](/internal-modes/classes-v2/discrete-transforms/improjection/columncount.html) Number of output coefficient coordinates.
  + [`columnLabels`](/internal-modes/classes-v2/discrete-transforms/improjection/columnlabels.html) Scientific column labels, distinct from array indices.
  + [`forwardMatrix`](/internal-modes/classes-v2/discrete-transforms/improjection/forwardmatrix.html) Matrix mapping sample columns to basis coefficients.
  + [`gramMatrix`](/internal-modes/classes-v2/discrete-transforms/improjection/grammatrix.html) Sampled Gram or explicitly prescribed signed coefficient system.
  + [`majorantGramMatrix`](/internal-modes/classes-v2/discrete-transforms/improjection/majorantgrammatrix.html) Positive coefficient metric used for error magnitudes.
  + [`metricMatrix`](/internal-modes/classes-v2/discrete-transforms/improjection/metricmatrix.html) Signed observation metric, or empty for a prescribed dual.
  + [`projectionKind`](/internal-modes/classes-v2/discrete-transforms/improjection/projectionkind.html) Galerkin construction or a caller-prescribed physical source dual.
  + [`provenance`](/internal-modes/classes-v2/discrete-transforms/improjection/provenance.html) Caller-supplied construction identity and numerical provenance.
  + [`sampleCount`](/internal-modes/classes-v2/discrete-transforms/improjection/samplecount.html) Number of sampled source coordinates accepted by project.
  + [`samplePairingMatrix`](/internal-modes/classes-v2/discrete-transforms/improjection/samplepairingmatrix.html) Signed map from source sample values to coefficient pairings.
  + [`sampledBasis`](/internal-modes/classes-v2/discrete-transforms/improjection/sampledbasis.html) Sampled basis, or empty for a prescribed dual without synthesis.
  + [`supportsGramAssessment`](/internal-modes/classes-v2/discrete-transforms/improjection/supportsgramassessment.html) Whether an actual sampled basis defines a Gram comparison.
  + [`supportsRoundTripAssessment`](/internal-modes/classes-v2/discrete-transforms/improjection/supportsroundtripassessment.html) Whether synthesis data defines an active-column round trip.
  + [`targetGramIsPositiveDefinite`](/internal-modes/classes-v2/discrete-transforms/improjection/targetgramispositivedefinite.html) Whether the continuous target is positive definite on active columns.
  + [`targetGramMatrix`](/internal-modes/classes-v2/discrete-transforms/improjection/targetgrammatrix.html) Continuous signed Gram matrix in the supplied column coordinates.
+ Measure projection quality
  + [`gramConditionNumber`](/internal-modes/classes-v2/discrete-transforms/improjection/gramconditionnumber.html) Condition number of the active sampled Gram matrix.
  + [`gramError`](/internal-modes/classes-v2/discrete-transforms/improjection/gramerror.html) Magnitude-normalized active Gram operator discrepancy.
  + [`inverseMatrixConditionNumber`](/internal-modes/classes-v2/discrete-transforms/improjection/inversematrixconditionnumber.html) Condition number of the active sampled basis, or NaN without synthesis.
  + [`leakage`](/internal-modes/classes-v2/discrete-transforms/improjection/leakage.html) Measure projected check profiles relative to their positive norms.
  + [`productError`](/internal-modes/classes-v2/discrete-transforms/improjection/producterror.html) Measure product projection errors with a positive coefficient norm.
  + [`roundTripError`](/internal-modes/classes-v2/discrete-transforms/improjection/roundtriperror.html) Operator norm of the active round-trip error, or NaN without synthesis.
  + [`sampledGramRank`](/internal-modes/classes-v2/discrete-transforms/improjection/sampledgramrank.html) Numerical rank of the active sampled Gram matrix.
+ Apply projections
  + [`project`](/internal-modes/classes-v2/discrete-transforms/improjection/project.html) sampled profiles with the signed physical pairing.


---