---
layout: default
title: IMGeostrophicTransform
has_children: false
has_toc: false
mathjax: true
parent: Discrete transforms
grand_parent: Class documentation V2
nav_order: 7
---

#  IMGeostrophicTransform

Compose APV and zero-APV transforms at positive wavenumber.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMGeostrophicTransform</code></pre></div></div>

## Overview

For APV eigendepth $$h_j$$ and horizontal wavenumber $$\kappa$$,

$$
\mu_\kappa^j=\kappa^2+\frac{f_0^2}{g h_j}.
$$

The APV endpoint-response columns are

$$
\mathbf r_q^{\kappa j}=-\frac{f_0}{g\mu_\kappa^j}
\begin{bmatrix}B_{\mathrm s}^j\\G_j(z_b)\end{bmatrix},
$$

restricted to active endpoints, where
$$B_{\mathrm s}=G(0)-F(0)$$ for a free surface and
$$B_{\mathrm s}=G(0)$$ for a rigid lid.

Sampled APV and source inputs have leading dimensions
`nZ x nK x ...`. APV coefficient outputs have leading dimensions
`nAPVModes x nK x ...`, while endpoint and zero-APV arrays have
leading dimensions `nEndpoints x nK x ...`. Trailing field
dimensions and complex values are preserved.

```matlab
transform = IMGeostrophicTransform(apvTransform=apvTransform,zeroAPVModes=zeroAPVModes,g0=g0,gd=gd);
[Aq,A0] = transform.transformStateForward(APV=q,endpointAnomalies=b);
```




## Topics
+ Create geostrophic transforms
  + [`IMGeostrophicTransform`](/internal-modes/classes-v2/discrete-transforms/imgeostrophictransform/imgeostrophictransform.html) Create an APV and zero-APV composition transform.
+ Transform admissible states
  + [`transformStateBack`](/internal-modes/classes-v2/discrete-transforms/imgeostrophictransform/transformstateback.html) Reconstruct sampled APV and active endpoint anomalies.
  + [`transformStateForward`](/internal-modes/classes-v2/discrete-transforms/imgeostrophictransform/transformstateforward.html) Transform an admissible APV and endpoint-anomaly state.
+ Project generic sources
  + [`transformSourceForward`](/internal-modes/classes-v2/discrete-transforms/imgeostrophictransform/transformsourceforward.html) Project generic vorticity and displacement sources.
+ Inspect geostrophic transforms
  + [`activeEndpoints`](/internal-modes/classes-v2/discrete-transforms/imgeostrophictransform/activeendpoints.html) Active endpoint coordinates in canonical order.
  + [`apvEndpointResponse`](/internal-modes/classes-v2/discrete-transforms/imgeostrophictransform/apvendpointresponse.html) APV endpoint-response pages $$\mathsf R_q^\kappa$$.
  + [`compatibilityDiagnostics`](/internal-modes/classes-v2/discrete-transforms/imgeostrophictransform/compatibilitydiagnostics.html) Compatibility and singularity diagnostics.
  + [`f0`](/internal-modes/classes-v2/discrete-transforms/imgeostrophictransform/f0.html) Coriolis parameter $$f_0$$.
  + [`g0`](/internal-modes/classes-v2/discrete-transforms/imgeostrophictransform/g0.html) Surface generalized-energy acceleration $$g_0$$.
  + [`gd`](/internal-modes/classes-v2/discrete-transforms/imgeostrophictransform/gd.html) Bottom generalized-energy acceleration $$g_d$$.
  + [`k`](/internal-modes/classes-v2/discrete-transforms/imgeostrophictransform/k.html) Horizontal-wavenumber pages $$\kappa$$.


---