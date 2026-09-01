---
layout: default
title: IMInternalModesDiscreteTransformAssessment
has_children: false
has_toc: false
mathjax: true
parent: Discrete transforms
grand_parent: Class documentation V2
nav_order: 5
---

#  IMInternalModesDiscreteTransformAssessment

Store retained-band diagnostics for aligned F/G transforms.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMInternalModesDiscreteTransformAssessment</code></pre></div></div>

## Overview

One family-specific point-and-weight rule is assessed for every
leading family prefix. Certified construction records independently
refitted candidate counts in `certificationSearch`; its final
`weightFitModeCount` equals `retainedModeCount`. `prefixDiagnostics`
reports the worst requested channel for
each policy, while `variablePrefixDiagnostics` exposes the underlying
per-variable Gram, rank, round-trip, and conditioning records.
The combined table includes `modeCount`, `lastModeNumber`, the Gram
error and limiting variable, leakage error and limiting rejected mode,
coupled quadratic error and limiting product channel/source labels,
and cumulative per-policy and combined acceptance flags.
Policy structs expose `tolerance`, per-prefix `error` and `accepted`
arrays, `maximumAcceptedModeCount`, `limitingValue`, and `reason`.
Internal-mode quadratic products are projected with the signed
Pontryagin pairing and signed Gram solve. Their relative error is
measured with the induced positive Hilbert majorant, recorded as
`quadraticAliasingPolicy.errorNorm`. Signed rejected-mode leakage
and MDA coupled quadratic products remain unsupported.




## Topics
+ Other
  + [`actualPointCount`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransformassessment/actualpointcount.html) Actual physical point count.
  + [`availableVariables`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransformassessment/availablevariables.html) Direct forward channels assessed in canonical order.
  + [`candidateModeCount`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransformassessment/candidatemodecount.html) Number of candidate family columns.
  + [`candidateModeNumber`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransformassessment/candidatemodenumber.html) Physical labels for candidate family columns.
  + [`candidateTransform`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransformassessment/candidatetransform.html) Full candidate-family transform on the fixed rule.
  + [`certificationSearch`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransformassessment/certificationsearch.html) Independently refitted count-selection attempts.
  + [`gramPolicy`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransformassessment/grampolicy.html) Worst-channel normalized-Gram policy result.
  + [`gridDesign`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransformassessment/griddesign.html) Provenance of the physical sample grid.
  + [`leakagePolicy`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransformassessment/leakagepolicy.html) Worst-channel rejected-mode leakage policy result.
  + [`limitingPolicy`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransformassessment/limitingpolicy.html) Policy responsible for the retained count.
  + [`limitingVariable`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransformassessment/limitingvariable.html) Variable responsible for the combined retained limit.
  + [`prefixDiagnostics`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransformassessment/prefixdiagnostics.html) Worst-channel and coupled-product prefix diagnostics.
  + [`quadraticAliasingPolicy`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransformassessment/quadraticaliasingpolicy.html) Coupled quadratic policy: signed projection, majorant error norm.
  + [`requestedPointCount`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransformassessment/requestedpointcount.html) Exact requested physical point count.
  + [`retainedModeCount`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransformassessment/retainedmodecount.html) Number of commonly retained family columns.
  + [`retainedModeNumber`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransformassessment/retainedmodenumber.html) Physical labels for retained family columns.
  + [`retentionReason`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransformassessment/retentionreason.html) Readable retained-band explanation.
  + [`transform`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransformassessment/transform.html) Production transform for the common accepted family prefix.
  + [`variablePrefixDiagnostics`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransformassessment/variableprefixdiagnostics.html) Return the per-prefix diagnostic table for F or G.
  + [`weightFit`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransformassessment/weightfit.html) Shared weight fit, or empty for caller-supplied weights.
  + [`weightFitModeCount`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransformassessment/weightfitmodecount.html) Number of family modes used to fit the stored weights.
  + [`withCertificationMetadata`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransformassessment/withcertificationmetadata.html) Return this assessment with grid and count-search provenance.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + Construct assessments
    + [`IMInternalModesDiscreteTransformAssessment`](/internal-modes/classes-v2/discrete-transforms/iminternalmodesdiscretetransformassessment/iminternalmodesdiscretetransformassessment.html) Create an aligned internal-mode transform assessment.


---