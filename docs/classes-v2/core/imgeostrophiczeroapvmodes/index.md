---
layout: default
title: IMGeostrophicZeroAPVModes
has_children: false
has_toc: false
mathjax: true
parent: Core
grand_parent: Class documentation V2
nav_order: 9
---

#  IMGeostrophicZeroAPVModes

Describe canonical geostrophic zero-APV modes at fixed wavenumber.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMGeostrophicZeroAPVModes</code></pre></div></div>

## Overview

`IMGeostrophicZeroAPVModes` stores the coefficient-independent
boundary-value problem

$$
\frac{\partial}{\partial z}\left(\frac{f_0^2}{N^2(z)}\frac{\partial F}{\partial z}\right)-k^2F=0,
\qquad
G(z)=-\frac{g}{N^2(z)}\frac{\partial F}{\partial z}.
$$

The canonical response vector contains the requested subset of

$$
B_{\mathrm s}[F]=
\begin{cases}
G(0)-F(0), & \texttt{freeSurface},\\
G(0), & \texttt{rigidLid},
\end{cases}
\qquad
B_{\mathrm d}[F]=G(z_b).
$$

These responses are proportional to the physical endpoint anomaly
functionals:

$$
\eta_0[F]=\frac{f_0}{g}B_{\mathrm s}[F],
\qquad
\eta_d[F]=\frac{f_0}{g}B_{\mathrm d}[F].
$$
Equivalently, the physical endpoint functionals are

$$
\eta_0[F]=
\begin{cases}
-\dfrac{f_0}{N^2(0)}\dfrac{\partial F}{\partial z}(0)-\dfrac{f_0}{g}F(0), & \texttt{freeSurface},\\
-\dfrac{f_0}{N^2(0)}\dfrac{\partial F}{\partial z}(0), & \texttt{rigidLid},
\end{cases}
\qquad
\eta_d[F]= -\frac{f_0}{N^2(z_b)}\frac{\partial F}{\partial z}(z_b).
$$

Solvers return one mode with unit response at each requested endpoint
and zero response at the other endpoint. Generalized-energy
coefficients are supplied later to basis rotation methods; they do
not change this problem or decide which endpoint coordinates exist.
For two endpoints,

$$
\mathbf B[F_0^{\mathrm{sur}}]=(1,0)^T,
\qquad
\mathbf B[F_0^{\mathrm{bot}}]=(0,1)^T.
$$

```matlab
problem = IMGeostrophicZeroAPVModes.atWavenumber(N2=N2,zDomain=[-4000 0],f0=1e-4,k=1e-4,endpoints=["surface","bottom"],surfaceBoundary="freeSurface");
solver = IMSolverSpectral(nEVP=128);
boundaryModes = solver.solveGeostrophicZeroAPVModes(problem);
depthModes = boundaryModes.rotateBoundaryDepth(g0=-0.035,gd=0);
```




## Topics
+ Create geostrophic zero-APV problems
  + [`IMGeostrophicZeroAPVModes`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodes/imgeostrophiczeroapvmodes.html) Create a canonical geostrophic zero-APV problem.
  + [`atWavenumber`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodes/atwavenumber.html) Create geostrophic zero-APV modes at fixed wavenumber.
+ Summarize geostrophic zero-APV problems
  + [`summarize`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodes/summarize.html) Print a readable problem summary.
+ Inspect geostrophic zero-APV problems
  + [`N2`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodes/n2.html) Buoyancy frequency squared function.
  + [`endpoints`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodes/endpoints.html) Canonically ordered endpoint coordinates.
  + [`f0`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodes/f0.html) Coriolis parameter.
  + [`g`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodes/g.html) Gravitational acceleration.
  + [`k`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodes/k.html) Horizontal wavenumbers.
  + [`metadata`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodes/metadata.html) Additional creation metadata.
  + [`modesPerWavenumber`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodes/modesperwavenumber.html) Return the number of canonical modes for each wavenumber.
  + [`surfaceBoundary`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodes/surfaceboundary.html) Surface endpoint convention.
  + [`zDomain`](/internal-modes/classes-v2/core/imgeostrophiczeroapvmodes/zdomain.html) Physical vertical domain.


---