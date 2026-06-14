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

Store one scalar canonical endpoint condition.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMBoundaryCondition</code></pre></div></div>

## Overview

`IMBoundaryCondition` represents
$$-[a u-b(pu')]=\lambda[c u-d(pu')].$$
The coefficients are endpoint scalars. The EVP supplies the endpoint
location when it assembles matrix rows or computes metric weights.




## Topics
+ Create boundary conditions
  + [`IMBoundaryCondition`](/internal-modes/classes-v2/core/imboundarycondition/imboundarycondition.html) Create a scalar endpoint condition.
  + [`dirichlet`](/internal-modes/classes-v2/core/imboundarycondition/dirichlet.html) Create `u=0`.
  + [`neumann`](/internal-modes/classes-v2/core/imboundarycondition/neumann.html) Create `p*u_z=0`.
  + [`robin`](/internal-modes/classes-v2/core/imboundarycondition/robin.html) Create `a*u-b*p*u_z=0`.
+ Inspect boundary conditions
  + [`a`](/internal-modes/classes-v2/core/imboundarycondition/a.html) Coefficient multiplying endpoint value on the left.
  + [`b`](/internal-modes/classes-v2/core/imboundarycondition/b.html) Coefficient multiplying endpoint flux on the left.
  + [`c`](/internal-modes/classes-v2/core/imboundarycondition/c.html) Coefficient multiplying endpoint value on the eigenvalue side.
  + [`d`](/internal-modes/classes-v2/core/imboundarycondition/d.html) Coefficient multiplying endpoint flux on the eigenvalue side.
  + [`determinant`](/internal-modes/classes-v2/core/imboundarycondition/determinant.html) Return the signed endpoint determinant.
  + [`isEigenvalueDependent`](/internal-modes/classes-v2/core/imboundarycondition/iseigenvaluedependent.html) Return true when the eigenvalue side is active.
  + [`metricWeight`](/internal-modes/classes-v2/core/imboundarycondition/metricweight.html) Return the endpoint metric weight.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`endpointNumeratorMatrix`](/internal-modes/classes-v2/core/imboundarycondition/endpointnumeratormatrix.html) Return the active-endpoint numerator matrix.
  + [`robinEnergyCoefficient`](/internal-modes/classes-v2/core/imboundarycondition/robinenergycoefficient.html) Return the ordinary Robin endpoint quadratic coefficient.


---