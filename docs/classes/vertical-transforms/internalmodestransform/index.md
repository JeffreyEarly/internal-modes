---
layout: default
title: InternalModesTransform
has_children: false
has_toc: false
mathjax: true
parent: Vertical transforms
grand_parent: Class documentation
nav_order: 2
---

#  InternalModesTransform

Apply vertical modal transforms and compute vertical modal spectra.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef InternalModesTransform < CAAnnotatedClass</code></pre></div></div>

## Overview

`InternalModesTransform` stores forward and inverse vertical
operators for one `InternalModesBasis`. It works only on vertical
scalar fields or modal coefficients; full wave-vortex state
projectors remain outside this package.

Public matrices use the full transform grid `z`. Some components are
solved on component-specific active rows before being expanded back to
that grid. For example, rigid-lid G transforms omit the zero Dirichlet
endpoint rows during the solve, then store zero endpoint columns in the
forward matrix and zero endpoint rows in the inverse matrix.

For modal coefficients $$p_j$$ and $$q_j$$, the vertical
cross-spectrum is

$$
S[p,q]_j = m\,s_j\,\Re\{p_j q_j^*\},
$$

where $$s_j$$ is the component spectral weight and $$m$$ is an
optional horizontal multiplicity supplied by the caller.

```matlab
coefficients = transform.projectVertical(profile,component="G");
spectrum = transform.spectrum(coefficients,component="G");
```




## Topics
+ Inspect transform properties
  + [`D`](/internal-modes/classes/vertical-transforms/internalmodestransform/d.html) Water-column depth used by the barotropic F spectral weight.
  + [`N2`](/internal-modes/classes/vertical-transforms/internalmodestransform/n2.html) Buoyancy frequency squared sampled at `z`.
  + [`componentRoleF`](/internal-modes/classes/vertical-transforms/internalmodestransform/componentrolef.html) F component role label.
  + [`componentRoleG`](/internal-modes/classes/vertical-transforms/internalmodestransform/componentroleg.html) G component role label.
  + [`conditionNumberF`](/internal-modes/classes/vertical-transforms/internalmodestransform/conditionnumberf.html) Condition-number diagnostic for the F forward operator.
  + [`conditionNumberG`](/internal-modes/classes/vertical-transforms/internalmodestransform/conditionnumberg.html) Condition-number diagnostic for the G forward operator.
  + [`f0`](/internal-modes/classes/vertical-transforms/internalmodestransform/f0.html) Coriolis parameter used by the source mode solve.
  + [`forwardF`](/internal-modes/classes/vertical-transforms/internalmodestransform/forwardf.html) Forward F projection matrix from samples to modal coefficients.
  + [`forwardG`](/internal-modes/classes/vertical-transforms/internalmodestransform/forwardg.html) Forward G projection matrix from samples to modal coefficients.
  + [`forwardProjectionAvailableF`](/internal-modes/classes/vertical-transforms/internalmodestransform/forwardprojectionavailablef.html) True when a canonical F projection exists for this basis.
  + [`forwardProjectionAvailableG`](/internal-modes/classes/vertical-transforms/internalmodestransform/forwardprojectionavailableg.html) True when a canonical G projection exists for this basis.
  + [`g`](/internal-modes/classes/vertical-transforms/internalmodestransform/g.html) Gravitational acceleration used by the source mode solve.
  + [`gramErrorF`](/internal-modes/classes/vertical-transforms/internalmodestransform/gramerrorf.html) Relative round-trip error for the retained F modes.
  + [`gramErrorG`](/internal-modes/classes/vertical-transforms/internalmodestransform/gramerrorg.html) Relative round-trip error for the retained G modes.
  + [`hF`](/internal-modes/classes/vertical-transforms/internalmodestransform/hf.html) Equivalent depths for retained F modes.
  + [`hG`](/internal-modes/classes/vertical-transforms/internalmodestransform/hg.html) Equivalent depths for retained G modes.
  + [`inverseF`](/internal-modes/classes/vertical-transforms/internalmodestransform/inversef.html) Inverse F reconstruction matrix from coefficients to samples.
  + [`inverseG`](/internal-modes/classes/vertical-transforms/internalmodestransform/inverseg.html) Inverse G reconstruction matrix from coefficients to samples.
  + [`kappa`](/internal-modes/classes/vertical-transforms/internalmodestransform/kappa.html) Horizontal wavenumber metadata.
  + [`nonlinearAliasLimit`](/internal-modes/classes/vertical-transforms/internalmodestransform/nonlinearaliaslimit.html) Number of modes allowed by the nonlinear aliasing policy.
  + [`omega`](/internal-modes/classes/vertical-transforms/internalmodestransform/omega.html) Frequency metadata.
  + [`problemType`](/internal-modes/classes/vertical-transforms/internalmodestransform/problemtype.html) Text label identifying the source vertical eigenvalue problem.
  + [`projectionResolvedModes`](/internal-modes/classes/vertical-transforms/internalmodestransform/projectionresolvedmodes.html) Number of modes accepted by the projection-quality diagnostic.
  + [`rejectedModesF`](/internal-modes/classes/vertical-transforms/internalmodestransform/rejectedmodesf.html) Rejected F mode numbers.
  + [`rejectedModesG`](/internal-modes/classes/vertical-transforms/internalmodestransform/rejectedmodesg.html) Rejected G mode numbers.
  + [`retainedModesF`](/internal-modes/classes/vertical-transforms/internalmodestransform/retainedmodesf.html) Retained F mode numbers.
  + [`retainedModesG`](/internal-modes/classes/vertical-transforms/internalmodestransform/retainedmodesg.html) Retained G mode numbers.
  + [`selectionReason`](/internal-modes/classes/vertical-transforms/internalmodestransform/selectionreason.html) Prefix-selection reason for model-style transforms.
  + [`sourceDescription`](/internal-modes/classes/vertical-transforms/internalmodestransform/sourcedescription.html) Text label identifying the source solver or factory path.
  + [`transformStatusF`](/internal-modes/classes/vertical-transforms/internalmodestransform/transformstatusf.html) F transform construction status.
  + [`transformStatusG`](/internal-modes/classes/vertical-transforms/internalmodestransform/transformstatusg.html) G transform construction status.
  + [`weightsF`](/internal-modes/classes/vertical-transforms/internalmodestransform/weightsf.html) Quadrature weights used to build F numerical projections.
  + [`weightsG`](/internal-modes/classes/vertical-transforms/internalmodestransform/weightsg.html) Quadrature weights used to build G numerical projections.
  + [`z`](/internal-modes/classes/vertical-transforms/internalmodestransform/z.html) Depth grid for vertical fields acted on by this transform.
+ Apply vertical transforms
  + [`crossTransformTo`](/internal-modes/classes/vertical-transforms/internalmodestransform/crosstransformto.html) Map coefficients from a target transform into this transform.
  + [`forward`](/internal-modes/classes/vertical-transforms/internalmodestransform/forward.html) Return a forward vertical projection matrix.
  + [`inverse`](/internal-modes/classes/vertical-transforms/internalmodestransform/inverse.html) Return an inverse vertical reconstruction matrix.
  + [`projectVertical`](/internal-modes/classes/vertical-transforms/internalmodestransform/projectvertical.html) Project vertical samples onto modal coefficients.
  + [`reconstructVertical`](/internal-modes/classes/vertical-transforms/internalmodestransform/reconstructvertical.html) Reconstruct vertical samples from modal coefficients.
+ Analyze vertical spectra
  + [`crossSpectrum`](/internal-modes/classes/vertical-transforms/internalmodestransform/crossspectrum.html) Compute a vertical modal cross-spectrum.
  + [`spectralWeightsF`](/internal-modes/classes/vertical-transforms/internalmodestransform/spectralweightsf.html) Parseval weights used for F modal spectra.
  + [`spectralWeightsG`](/internal-modes/classes/vertical-transforms/internalmodestransform/spectralweightsg.html) Parseval weights used for G modal spectra.
  + [`spectrum`](/internal-modes/classes/vertical-transforms/internalmodestransform/spectrum.html) Compute a vertical modal auto-spectrum.
+ Persist vertical transforms
  + [`InternalModesTransform`](/internal-modes/classes/vertical-transforms/internalmodestransform/internalmodestransform.html) Create a vertical transform from canonical persisted state.
+ Other
  + [`classRequiredPropertyNames`](/internal-modes/classes/vertical-transforms/internalmodestransform/classrequiredpropertynames.html)


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`modeF`](/internal-modes/classes/vertical-transforms/internalmodestransform/modef.html) F-mode coordinate used by annotated NetCDF persistence.
  + [`modeG`](/internal-modes/classes/vertical-transforms/internalmodestransform/modeg.html) G-mode coordinate used by annotated NetCDF persistence.
  + [`withDiagnostics`](/internal-modes/classes/vertical-transforms/internalmodestransform/withdiagnostics.html) Return a copy with updated selection diagnostics.
  + [`zIndex`](/internal-modes/classes/vertical-transforms/internalmodestransform/zindex.html) Row coordinate used by annotated NetCDF persistence.


---