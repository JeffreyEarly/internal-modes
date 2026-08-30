---
layout: default
title: IMGeostrophicZeroAPVModesBasis
has_children: false
has_toc: false
mathjax: true
parent: Core
grand_parent: Class documentation V2
nav_order: 8
---

#  IMGeostrophicZeroAPVModesBasis

Store canonical or rotated geostrophic zero-APV modes.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMGeostrophicZeroAPVModesBasis</code></pre></div></div>

## Overview

The canonical columns have unit endpoint responses. For every
wavenumber page, the selected response matrix is the identity. The
diagnostic structure is

$$
G(z)=-\frac{g}{N^2(z)}\frac{\partial F}{\partial z}(z).
$$

`F(z)` and `G(z)` have dimensions `nZ x nEndpoints x nK`. The four
public quadratic-form properties have dimensions
`nEndpoints x nEndpoints x nK`. Rotation methods apply the same
boundary-coordinate matrix to both physical variables.

If $$\mathbf b_{\mathrm s}$$ and $$\mathbf b_{\mathrm d}$$ are the
surface and bottom response rows and $$\mathbf F_{\mathrm s}$$ and
$$\mathbf F_{\mathrm d}$$ are the corresponding value rows, then
each canonical wavenumber page stores

$$
\mathsf R_B=\mathbf B^T\mathbf B,
\qquad
\mathsf H=\frac{f_0^2}{2gk^4}
\operatorname{sym}\left(-\mathbf F_{\mathrm s}^T\mathbf b_{\mathrm s}+\mathbf F_{\mathrm d}^T\mathbf b_{\mathrm d}\right),
$$

$$
\mathsf B_0=\frac{f_0^2}{2g^2k^4}\mathbf b_{\mathrm s}^T\mathbf b_{\mathrm s},
\qquad
\mathsf B_d=\frac{f_0^2}{2g^2k^4}\mathbf b_{\mathrm d}^T\mathbf b_{\mathrm d}.
$$

An absent endpoint has a zero form matrix. Generalized energy is
$$\mathsf H_g=\mathsf H+g_0\mathsf B_0+g_d\mathsf B_d$$.

```matlab
boundaryModes = solver.solveGeostrophicZeroAPVModes(problem);
depthModes = boundaryModes.rotateBoundaryDepth(g0=-0.035,gd=0);
surfaceModes = boundaryModes.rotateSurfaceBuoyancy(g0=-0.035,gd=0);
Hg = boundaryModes.generalizedEnergyMatrix(g0=-0.035,gd=0);
customModes = boundaryModes.rotateWithPencil(name="custom",leftMatrix=Hg,rightMatrix=boundaryModes.endpointResponseMetric);
```




## Topics
+ Evaluate geostrophic zero-APV modes
  + [`F`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/f.html) Evaluate streamfunction structures $$F(z)$$.
  + [`G`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/g.html) Evaluate diagnostic displacement structures $$G(z)$$.
  + [`IMGeostrophicZeroAPVModesBasis`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/imgeostrophiczeroapvmodesbasis.html) Create a canonical boundary-normalized zero-APV basis.
+ Inspect geostrophic zero-APV modes
  + [`N2`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/n2.html) Buoyancy frequency squared function.
  + [`endpoints`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/endpoints.html) Canonically ordered endpoint-coordinate labels.
  + [`f0`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/f0.html) Coriolis parameter.
  + [`g`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/g_.html) Gravitational acceleration.
  + [`h0`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/h0.html) Signed boundary depths.
  + [`k`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/k.html) Horizontal-wavenumber pages.
  + [`metadata`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/metadata.html) Additional creation and rotation metadata.
  + [`normalizationConvention`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/normalizationconvention.html) Rotation normalization convention.
  + [`problem`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/problem.html) Geostrophic zero-APV problem descriptor.
  + [`rotationEigenvalues`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/rotationeigenvalues.html) Pagewise generalized-pencil eigenvalues.
  + [`rotationMatrix`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/rotationmatrix.html) Canonical-to-current column transformation.
  + [`rotationName`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/rotationname.html) Name of the current boundary-coordinate rotation.
  + [`rotationResiduals`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/rotationresiduals.html) Pagewise pencil and normalization residuals.
  + [`signatures`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/signatures.html) Pagewise quadratic-form signatures.
  + [`solver`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/solver.html) Numerical solver used to compute the canonical structures.
  + [`summarize`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/summarize.html) Print a readable geostrophic zero-APV basis summary.
  + [`surfaceBoundary`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/surfaceboundary.html) Surface endpoint convention.
  + [`zDomain`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/zdomain.html) Physical vertical domain.
+ Form geostrophic zero-APV quadratic forms
  + [`bottomBuoyancyMatrix`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/bottombuoyancymatrix.html) Bottom-buoyancy matrix $$\mathsf B_d$$.
  + [`endpointResponseMetric`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/endpointresponsemetric.html) Endpoint-response Gram matrix $$\mathsf R_B$$.
  + [`energyMatrix`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/energymatrix.html) Physical-energy matrix $$\mathsf H$$.
  + [`generalizedEnergyMatrix`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/generalizedenergymatrix.html) Return $$\mathsf H_g=\mathsf H+g_0\mathsf B_0+g_d\mathsf B_d$$.
  + [`surfaceBuoyancyMatrix`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/surfacebuoyancymatrix.html) Surface-buoyancy matrix $$\mathsf B_0$$.
+ Rotate geostrophic zero-APV modes
  + [`rotateBoundaryDepth`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/rotateboundarydepth.html) Diagonalize generalized energy relative to endpoint response.
  + [`rotateSurfaceBuoyancy`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/rotatesurfacebuoyancy.html) Diagonalize surface buoyancy relative to generalized energy.
  + [`rotateWithPencil`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodesbasis/rotatewithpencil.html) Apply a custom symmetric matrix-pencil rotation.


---