---
layout: default
title: IMAnalyticalGeostrophicZeroAPVModesBasis
has_children: false
has_toc: false
mathjax: true
parent: Analytical bases
grand_parent: Class documentation V2
nav_order: 5
---

#  IMAnalyticalGeostrophicZeroAPVModesBasis

Store exact canonical or rotated geostrophic zero-APV modes.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMAnalyticalGeostrophicZeroAPVModesBasis</code></pre></div></div>

## Overview

Concrete analytical stratification families create the canonical
columns. Their endpoint responses satisfy

$$
\mathbf B[F_0^{\mathrm{sur}}]=(1,0)^T,
\qquad
\mathbf B[F_0^{\mathrm{bot}}]=(0,1)^T,
$$

for the requested subset of surface and bottom endpoints. The exact
diagnostic variable is

$$
G(z)=-\frac{g}{N^2(z)}\frac{\partial F}{\partial z}(z).
$$

`F(z)` and `G(z)` have dimensions `nZ x nEndpoints x nK`. Quadratic
forms have dimensions `nEndpoints x nEndpoints x nK`. Rotations apply
the same pagewise column map to both exact variables.

```matlab
solution = IMConstantStratificationSolution(N0=5.2e-3,zDomain=[-4000 0],f0=1e-4);
exactModes = solution.geostrophicZeroAPVModesAtWavenumber(1e-4);
depthModes = exactModes.rotateBoundaryDepth(g0=-0.03,gd=0.01);
```




## Topics
+ Evaluate exact geostrophic zero-APV modes
  + [`F`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/f.html) Evaluate exact streamfunction structures $$F(z)$$.
  + [`G`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/g.html) Evaluate exact diagnostic displacement structures $$G(z)$$.
+ Inspect exact geostrophic zero-APV modes
  + [`N2`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/n2.html) Buoyancy frequency squared function.
  + [`endpoints`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/endpoints.html) Canonically ordered endpoint-coordinate labels.
  + [`f0`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/f0.html) Coriolis parameter.
  + [`g`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/g_.html) Gravitational acceleration.
  + [`h0`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/h0.html) Signed boundary depths from the boundary-depth rotation.
  + [`k`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/k.html) Horizontal-wavenumber pages.
  + [`metadata`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/metadata.html) Additional creation and rotation metadata.
  + [`normalizationConvention`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/normalizationconvention.html) Rotation normalization convention.
  + [`problem`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/problem.html) Canonical geostrophic zero-APV problem represented by the formulas.
  + [`rotationEigenvalues`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/rotationeigenvalues.html) Pagewise generalized-pencil eigenvalues.
  + [`rotationMatrix`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/rotationmatrix.html) Canonical-to-current pagewise column transformation.
  + [`rotationName`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/rotationname.html) Name of the current boundary-coordinate rotation.
  + [`rotationResiduals`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/rotationresiduals.html) Pagewise pencil and normalization residuals.
  + [`signatures`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/signatures.html) Pagewise quadratic-form signatures.
  + [`solution`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/solution.html) Analytical solution family that created this basis.
  + [`summarize`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/summarize.html) Print a readable exact zero-APV basis summary.
  + [`surfaceBoundary`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/surfaceboundary.html) Surface endpoint convention.
  + [`zDomain`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/zdomain.html) Physical vertical domain.
+ Form exact geostrophic zero-APV quadratic forms
  + [`bottomBuoyancyMatrix`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/bottombuoyancymatrix.html) Bottom-buoyancy matrix $$\mathsf B_d$$.
  + [`endpointResponseMetric`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/endpointresponsemetric.html) Endpoint-response Gram matrix $$\mathsf R_B$$.
  + [`energyMatrix`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/energymatrix.html) Physical-energy matrix $$\mathsf H$$.
  + [`generalizedEnergyMatrix`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/generalizedenergymatrix.html) Return $$\mathsf H_g=\mathsf H+g_0\mathsf B_0+g_d\mathsf B_d$$.
  + [`surfaceBuoyancyMatrix`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/surfacebuoyancymatrix.html) Surface-buoyancy matrix $$\mathsf B_0$$.
+ Rotate exact geostrophic zero-APV modes
  + [`rotateBoundaryDepth`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/rotateboundarydepth.html) Diagonalize generalized energy relative to endpoint response.
  + [`rotateSurfaceBuoyancy`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/rotatesurfacebuoyancy.html) Diagonalize surface buoyancy relative to generalized energy.
  + [`rotateWithPencil`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/rotatewithpencil.html) Apply a custom symmetric matrix-pencil rotation.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`IMAnalyticalGeostrophicZeroAPVModesBasis`](/internal-modes/classes-v2/analytical-bases/imanalyticalgeostrophiczeroapvmodesbasis/imanalyticalgeostrophiczeroapvmodesbasis.html) Create an exact canonical boundary-normalized basis.


---