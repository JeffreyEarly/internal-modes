---
layout: default
title: IMBoundary
has_children: false
has_toc: false
mathjax: true
parent: Core
grand_parent: Class documentation V2
nav_order: 4
---

#  IMBoundary

Describe an internal-mode boundary law.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMBoundary</code></pre></div></div>

## Overview

`IMBoundary` stores the physical boundary law supplied to an
`IMEigenvalueProblem`. Standard EVPs keep one law at the surface and
one law at the bottom. When an EVP is assembled, each law is resolved
with the EVP formulation and endpoint to produce the operator row,
endpoint weights, and mode-index metadata needed by the solver and
basis-set layers.

Boundary operators use the physical coordinate derivative
$$\partial_z$$. A boundary functional has the form
$$\sum_p a_p(z_\ell)\partial_z^p u(z_\ell)$$ at endpoint
$$z_\ell$$; no outward-normal derivative sign is inserted at the
bottom.

Boundary weights are not boundary laws. They are endpoint products
added to modal inner products so the EVP remains orthogonal under the
chosen boundary law. Endpoint orientation enters through Green's
identity:
$$s_\ell=+1$$ at the surface and $$s_\ell=-1$$ at the bottom, so a
generated term with coefficient $$c$$ contributes
$$s_\ell c\,u_i(z_\ell)v_j(z_\ell)$$.

```matlab
evp = IMEigenvalueProblem.hydrostaticGModes( ...
    surfaceBoundary=IMBoundary.free(), bottomBoundary=IMBoundary.rigid());
```




## Topics
+ Create boundary laws
  + [`custom`](/internal-modes/classes-v2/core/imboundary/custom.html) Create a custom location-free operator boundary law.
  + [`dirichlet`](/internal-modes/classes-v2/core/imboundary/dirichlet.html) Create a location-free homogeneous Dirichlet boundary law.
  + [`free`](/internal-modes/classes-v2/core/imboundary/free.html) Create a location-free free boundary law.
  + [`linearF`](/internal-modes/classes-v2/core/imboundary/linearf.html) Create a location-free linear `F` boundary law.
  + [`linearG`](/internal-modes/classes-v2/core/imboundary/linearg.html) Create a location-free linear `G` boundary law.
  + [`neumann`](/internal-modes/classes-v2/core/imboundary/neumann.html) Create a location-free homogeneous Neumann boundary law.
  + [`noSlip`](/internal-modes/classes-v2/core/imboundary/noslip.html) Create a location-free no-slip boundary law.
  + [`rigid`](/internal-modes/classes-v2/core/imboundary/rigid.html) Create a location-free rigid boundary law.
+ Inspect boundary laws
  + [`boundaryModeNumber`](/internal-modes/classes-v2/core/imboundary/boundarymodenumber.html) Physical mode number for an endpoint boundary branch.
  + [`boundaryWeights`](/internal-modes/classes-v2/core/imboundary/boundaryweights.html) Endpoint weights implied by this law for modal inner products.
  + [`family`](/internal-modes/classes-v2/core/imboundary/family.html) Boundary family identifier.
  + [`hasKnownBoundaryWeights`](/internal-modes/classes-v2/core/imboundary/hasknownboundaryweights.html) True when the compatible boundary weights are known.
  + [`indexRank`](/internal-modes/classes-v2/core/imboundary/indexrank.html) Rank of the boundary-mode eigenvalue contribution.
  + [`indexSign`](/internal-modes/classes-v2/core/imboundary/indexsign.html) Sign of the boundary-mode eigenvalue contribution.
  + [`leftOperator`](/internal-modes/classes-v2/core/imboundary/leftoperator.html) Left-side boundary functional.
  + [`location`](/internal-modes/classes-v2/core/imboundary/location.html) Physical endpoint location.
  + [`rightOperator`](/internal-modes/classes-v2/core/imboundary/rightoperator.html) Right-side eigenvalue boundary functional.
  + [`variable`](/internal-modes/classes-v2/core/imboundary/variable.html) constrained by this boundary law.


---