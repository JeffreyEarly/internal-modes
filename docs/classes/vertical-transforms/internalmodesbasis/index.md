---
layout: default
title: InternalModesBasis
has_children: false
has_toc: false
mathjax: true
parent: Vertical transforms
grand_parent: Class documentation
nav_order: 1
---

#  InternalModesBasis

Store solved vertical modes together with their component roles.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef InternalModesBasis < CAAnnotatedClass</code></pre></div></div>

## Overview

`InternalModesBasis` is the canonical container for a solved vertical
mode basis. It records the sampled inverse modes $$F_j(z)$$ and
$$G_j(z)$$, the equivalent depths $$h_j$$, the quadrature grid, and
whether each component has a Sturm-Liouville forward projection.

Geostrophic bases at $$\omega=0$$ have canonical F and G projections,

$$
\mathcal{F}_g^j[u] = \gamma_j^{-1}\int_{-D}^{0} u F_g^j\,dz,
\qquad
\mathcal{G}_g^j[\eta] = \frac{1}{g}\int_{-D}^{0} N^2 \eta G_g^j\,dz.
$$

Nonzero-$$\kappa$$ IGW bases have a canonical G projection,

$$
\mathcal{G}_\kappa^j[\eta] =
\frac{1}{g}\int_{-D}^{0} \left(N^2-f_0^2\right)\eta G_\kappa^j\,dz,
$$

but no independent canonical wave-F projection.

```matlab
basis = InternalModesBasis.fromSolverAtFrequency(im,0,nModes=32);
transform = basis.nativeTransform(component="G");
```




## Topics
+ Create vertical bases
  + [`InternalModesBasis`](/internal-modes/classes/vertical-transforms/internalmodesbasis/internalmodesbasis.html) Create a persisted vertical basis from canonical state.
  + [`fromSolvedModes`](/internal-modes/classes/vertical-transforms/internalmodesbasis/fromsolvedmodes.html) Create a basis from already-solved mode arrays.
  + [`fromSolverAtFrequency`](/internal-modes/classes/vertical-transforms/internalmodesbasis/fromsolveratfrequency.html) Solve modes at fixed frequency and return an annotated basis.
  + [`fromSolverAtWavenumber`](/internal-modes/classes/vertical-transforms/internalmodesbasis/fromsolveratwavenumber.html) Solve modes at fixed horizontal wavenumber.
+ Inspect basis properties
  + [`D`](/internal-modes/classes/vertical-transforms/internalmodesbasis/d.html) Water-column depth $$D=z_{\max}-z_{\min}$$.
  + [`F`](/internal-modes/classes/vertical-transforms/internalmodesbasis/f.html) Sampled F inverse modes with rows matching `z`.
  + [`G`](/internal-modes/classes/vertical-transforms/internalmodesbasis/g.html) Sampled G inverse modes with rows matching `z`.
  + [`N2`](/internal-modes/classes/vertical-transforms/internalmodesbasis/n2.html) Buoyancy frequency squared sampled at `z`.
  + [`componentRoleF`](/internal-modes/classes/vertical-transforms/internalmodesbasis/componentrolef.html) Component role label for the F component.
  + [`componentRoleG`](/internal-modes/classes/vertical-transforms/internalmodesbasis/componentroleg.html) Component role label for the G component.
  + [`f0`](/internal-modes/classes/vertical-transforms/internalmodesbasis/f0.html) Coriolis parameter used by the mode solve.
  + [`forwardProjectionAvailableF`](/internal-modes/classes/vertical-transforms/internalmodesbasis/forwardprojectionavailablef.html) Boolean flag indicating whether canonical F projection exists.
  + [`forwardProjectionAvailableG`](/internal-modes/classes/vertical-transforms/internalmodesbasis/forwardprojectionavailableg.html) Boolean flag indicating whether canonical G projection exists.
  + [`g`](/internal-modes/classes/vertical-transforms/internalmodesbasis/g_.html) Gravitational acceleration used by the mode solve.
  + [`h`](/internal-modes/classes/vertical-transforms/internalmodesbasis/h.html) Equivalent depths associated with the stored mode columns.
  + [`kappa`](/internal-modes/classes/vertical-transforms/internalmodesbasis/kappa.html) Horizontal wavenumber associated with the solve.
  + [`omega`](/internal-modes/classes/vertical-transforms/internalmodesbasis/omega.html) Frequency associated with the solve.
  + [`orthogonalityWeightF`](/internal-modes/classes/vertical-transforms/internalmodesbasis/orthogonalityweightf.html) Orthogonality weight label for F projections.
  + [`orthogonalityWeightG`](/internal-modes/classes/vertical-transforms/internalmodesbasis/orthogonalityweightg.html) Orthogonality weight label for G projections.
  + [`problemType`](/internal-modes/classes/vertical-transforms/internalmodesbasis/problemtype.html) Text label identifying the vertical eigenvalue problem.
  + [`sourceDescription`](/internal-modes/classes/vertical-transforms/internalmodesbasis/sourcedescription.html) Text label identifying the source solver or factory path.
  + [`z`](/internal-modes/classes/vertical-transforms/internalmodesbasis/z.html) Depth grid where the inverse modes are sampled.
+ Build vertical transforms
  + [`fixedGridTransform`](/internal-modes/classes/vertical-transforms/internalmodesbasis/fixedgridtransform.html) Build a transform on another transform's vertical grid.
  + [`modelTransform`](/internal-modes/classes/vertical-transforms/internalmodesbasis/modeltransform.html) Build a prefix-retained transform for model vertical grids.
  + [`nativeTransform`](/internal-modes/classes/vertical-transforms/internalmodesbasis/nativetransform.html) Build a transform on the basis' native vertical grid.
  + [`observationProjection`](/internal-modes/classes/vertical-transforms/internalmodesbasis/observationprojection.html) Build a vertical projection for an arbitrary observation grid.
+ Analyze vertical spectra
  + [`spectralWeights`](/internal-modes/classes/vertical-transforms/internalmodesbasis/spectralweights.html) Return Parseval weights for modal spectra.
+ Other
  + [`classRequiredPropertyNames`](/internal-modes/classes/vertical-transforms/internalmodesbasis/classrequiredpropertynames.html)


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`modeF`](/internal-modes/classes/vertical-transforms/internalmodesbasis/modef.html) F-mode index coordinate used by annotated NetCDF persistence.
  + [`modeG`](/internal-modes/classes/vertical-transforms/internalmodesbasis/modeg.html) G-mode index coordinate used by annotated NetCDF persistence.
  + [`zIndex`](/internal-modes/classes/vertical-transforms/internalmodesbasis/zindex.html) Row coordinate used by annotated NetCDF persistence.


---