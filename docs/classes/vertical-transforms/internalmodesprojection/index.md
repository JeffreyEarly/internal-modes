---
layout: default
title: InternalModesProjection
has_children: false
has_toc: false
mathjax: true
parent: Vertical transforms
grand_parent: Class documentation
nav_order: 3
---

#  InternalModesProjection

Project arbitrary vertical observations onto resolvable modes.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef InternalModesProjection < CAAnnotatedClass</code></pre></div></div>

## Overview

`InternalModesProjection` stores the weighted least-squares operator
for an observation matrix $$B=H\Phi$$. Unlike model transforms, an
observation grid may identify a non-contiguous subset of modes. The
resolution matrix

$$
R = A B
$$

describes how true candidate-mode coefficients appear in the
recovered retained coefficients.

```matlab
projection = basis.observationProjection(zObs,component="G");
etaHat = projection.project(etaObs);
```




## Topics
+ Inspect projection properties
  + [`aliasingMatrix`](/internal-modes/classes/vertical-transforms/internalmodesprojection/aliasingmatrix.html) Alias matrix formed from rejected-mode columns of the resolution matrix.
  + [`candidateModes`](/internal-modes/classes/vertical-transforms/internalmodesprojection/candidatemodes.html) Candidate mode numbers considered by the projection.
  + [`component`](/internal-modes/classes/vertical-transforms/internalmodesprojection/component.html) label, either `"F"` or `"G"`.
  + [`componentRole`](/internal-modes/classes/vertical-transforms/internalmodesprojection/componentrole.html) Component role label.
  + [`conditionNumber`](/internal-modes/classes/vertical-transforms/internalmodesprojection/conditionnumber.html) Condition number of `gramMatrix`.
  + [`diagnosticFinalRelativePivot`](/internal-modes/classes/vertical-transforms/internalmodesprojection/diagnosticfinalrelativepivot.html) Final retained relative QR pivot.
  + [`diagnosticRankTolerance`](/internal-modes/classes/vertical-transforms/internalmodesprojection/diagnosticranktolerance.html) QR rank tolerance used for mode selection.
  + [`gramMatrix`](/internal-modes/classes/vertical-transforms/internalmodesprojection/grammatrix.html) Retained Gram matrix used in the weighted least-squares solve.
  + [`isCanonicalComponent`](/internal-modes/classes/vertical-transforms/internalmodesprojection/iscanonicalcomponent.html) True when this projection has a canonical modal spectrum.
  + [`observationMatrix`](/internal-modes/classes/vertical-transforms/internalmodesprojection/observationmatrix.html) Explicit observation matrix, when supplied by the caller.
  + [`observationZ`](/internal-modes/classes/vertical-transforms/internalmodesprojection/observationz.html) Observation depths for point-sampled projections.
  + [`projectionMatrix`](/internal-modes/classes/vertical-transforms/internalmodesprojection/projectionmatrix.html) Forward projection matrix from observations to retained coefficients.
  + [`reconstructionMatrix`](/internal-modes/classes/vertical-transforms/internalmodesprojection/reconstructionmatrix.html) Reconstruction matrix from retained coefficients to native samples.
  + [`rejectedModes`](/internal-modes/classes/vertical-transforms/internalmodesprojection/rejectedmodes.html) Rejected candidate mode numbers.
  + [`resolutionMatrix`](/internal-modes/classes/vertical-transforms/internalmodesprojection/resolutionmatrix.html) Resolution matrix from true candidate coefficients to recovered rows.
  + [`retainedModes`](/internal-modes/classes/vertical-transforms/internalmodesprojection/retainedmodes.html) Retained mode numbers selected as resolvable.
  + [`transformStatus`](/internal-modes/classes/vertical-transforms/internalmodesprojection/transformstatus.html) Transform construction status.
+ Apply observation projections
  + [`project`](/internal-modes/classes/vertical-transforms/internalmodesprojection/project.html) observations onto retained modal coefficients.
  + [`reconstruct`](/internal-modes/classes/vertical-transforms/internalmodesprojection/reconstruct.html) sampled values from retained coefficients.
+ Analyze observation spectra
  + [`crossSpectrum`](/internal-modes/classes/vertical-transforms/internalmodesprojection/crossspectrum.html) Compute a retained-mode observation cross-spectrum.
  + [`expectedRecoveredSpectrum`](/internal-modes/classes/vertical-transforms/internalmodesprojection/expectedrecoveredspectrum.html) Apply the spectral window to a true candidate-mode spectrum.
  + [`spectralWeights`](/internal-modes/classes/vertical-transforms/internalmodesprojection/spectralweights.html) Candidate-mode spectral weights.
  + [`spectralWindow`](/internal-modes/classes/vertical-transforms/internalmodesprojection/spectralwindow.html) Squared resolution matrix used as an expected spectral window.
  + [`spectrum`](/internal-modes/classes/vertical-transforms/internalmodesprojection/spectrum.html) Compute a retained-mode observation spectrum.
+ Persist observation projections
  + [`InternalModesProjection`](/internal-modes/classes/vertical-transforms/internalmodesprojection/internalmodesprojection.html) Create an observation projection from canonical persisted state.
+ Other
  + [`classRequiredPropertyNames`](/internal-modes/classes/vertical-transforms/internalmodesprojection/classrequiredpropertynames.html)


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`candidateModeIndex`](/internal-modes/classes/vertical-transforms/internalmodesprojection/candidatemodeindex.html) Candidate-mode coordinate used by annotated persistence.
  + [`observationIndex`](/internal-modes/classes/vertical-transforms/internalmodesprojection/observationindex.html) Observation index coordinate used by annotated persistence.
  + [`retainedModeIndex`](/internal-modes/classes/vertical-transforms/internalmodesprojection/retainedmodeindex.html) Retained-mode coordinate used by annotated persistence.
  + [`selectResolvableModes`](/internal-modes/classes/vertical-transforms/internalmodesprojection/selectresolvablemodes.html) Select observation-resolvable columns with pivoted QR.
  + [`weightedProjection`](/internal-modes/classes/vertical-transforms/internalmodesprojection/weightedprojection.html) Build a weighted projection for selected observation columns.


---