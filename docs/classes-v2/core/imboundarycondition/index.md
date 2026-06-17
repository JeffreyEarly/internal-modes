---
layout: default
title: IMBoundaryCondition
has_children: false
has_toc: false
mathjax: true
parent: Core
grand_parent: Class documentation V2
nav_order: 2
---

#  IMBoundaryCondition

Store one scalar canonical boundary condition.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMBoundaryCondition</code></pre></div></div>

## Overview

`IMBoundaryCondition` represents
$$-[a u-b(pu')]=\lambda[c u-d(pu')].$$
The coefficients define the boundary condition. The EVP supplies the
endpoint location when it assembles matrix rows or computes endpoint weight
coefficients. The stored properties `a`, `b`, `c`, and `d` define
the boundary condition; `determinant(location)` derives the signed
determinant $$D_i$$; and `endpointWeightCoefficient(location)`
derives the scalar $$D_i^{-1}$$ used by
`IMEigenvalueProblem.endpointWeights`.

At endpoint $$z_\ell$$, define
$$P_\ell=p(z_\ell)\frac{\partial u}{\partial z}(z_\ell).$$
The signed determinant is
$$D_\ell=\sigma_\ell(ad-bc),\qquad
\sigma_\mathrm{bottom}=+1,\quad \sigma_\mathrm{surface}=-1.$$
Only active boundary conditions, where `(c,d) ~= (0,0)`, produce
endpoint norm weights:

| Boundary condition | Endpoint equation | Endpoint norm contribution |
| --- | --- | --- |
| Dirichlet | $$u_\ell=0$$ | none |
| Neumann | $$P_\ell=0$$ | none |
| Robin | $$a u_\ell-bP_\ell=0$$ | none; `robinEnergyCoefficient` is an energy diagnostic |
| Active value condition | eigenvalue side uses $$c u_\ell$$ | $$D_\ell^{-1}(c u_\ell)^2$$ |
| Active flux condition | eigenvalue side uses $$-dP_\ell$$ | $$D_\ell^{-1}(-dP_\ell)^2$$ |
| General active condition | $$-[a u_\ell-bP_\ell]=\lambda[c u_\ell-dP_\ell]$$ | $$D_\ell^{-1}(c u_\ell-dP_\ell)^2$$ |
| Degenerate active condition | `(c,d) ~= (0,0)` and $$D_\ell=0$$ | unavailable |

For example, `IMBoundaryCondition(a=0,b=1,c=1,d=0)` has
$$D_\mathrm{surface}=+1$$ and $$D_\mathrm{bottom}=-1$$, so the
same boundary condition contributes $$+u_s^2$$ at the surface and
$$-u_b^2$$ at the bottom.




## Topics
+ Create boundary conditions
  + [`IMBoundaryCondition`](/internal-modes/classes-v2/core/imboundarycondition/imboundarycondition.html) Create a scalar boundary condition.
  + [`dirichlet`](/internal-modes/classes-v2/core/imboundarycondition/dirichlet.html) Create `u=0`.
  + [`neumann`](/internal-modes/classes-v2/core/imboundarycondition/neumann.html) Create `p*u_z=0`.
  + [`robin`](/internal-modes/classes-v2/core/imboundarycondition/robin.html) Create `a*u-b*p*u_z=0`.
+ Inspect boundary conditions
  + [`a`](/internal-modes/classes-v2/core/imboundarycondition/a.html) Coefficient multiplying endpoint value on the left.
  + [`b`](/internal-modes/classes-v2/core/imboundarycondition/b.html) Coefficient multiplying endpoint flux on the left.
  + [`c`](/internal-modes/classes-v2/core/imboundarycondition/c.html) Coefficient multiplying endpoint value on the eigenvalue side.
  + [`d`](/internal-modes/classes-v2/core/imboundarycondition/d.html) Coefficient multiplying endpoint flux on the eigenvalue side.
  + [`determinant`](/internal-modes/classes-v2/core/imboundarycondition/determinant.html) Return the signed endpoint determinant.
  + [`endpointWeightCoefficient`](/internal-modes/classes-v2/core/imboundarycondition/endpointweightcoefficient.html) Return the endpoint weight coefficient.
  + [`isEigenvalueDependent`](/internal-modes/classes-v2/core/imboundarycondition/iseigenvaluedependent.html) Return true when the eigenvalue side is active.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`endpointNumeratorMatrix`](/internal-modes/classes-v2/core/imboundarycondition/endpointnumeratormatrix.html) Return the active-endpoint numerator matrix.
  + [`robinEnergyCoefficient`](/internal-modes/classes-v2/core/imboundarycondition/robinenergycoefficient.html) Return the ordinary Robin endpoint quadratic coefficient.


---