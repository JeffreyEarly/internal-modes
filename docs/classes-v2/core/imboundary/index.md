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

Describe an internal-mode boundary condition.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMBoundary</code></pre></div></div>

## Overview

`IMBoundary` stores a mathematical boundary condition. A boundary can
be location-free, such as `IMBoundary.free()`, or placed at a physical
endpoint, such as `IMBoundary.free().at("surface",formulation="G")`.
EVP factories normally accept location-free surface and bottom laws
and place them through `IMBoundary.conditions(...)`.

Boundary operators use the physical coordinate derivative
$$\partial_z$$ at both endpoints. A boundary functional has the form
$$\sum_p a_p(z_\ell)\partial_z^p u(z_\ell)$$ at endpoint
$$z_\ell$$; no outward-normal derivative sign is inserted at the
bottom.

Boundary inner-product terms are not boundary conditions. They are the
boundary trace products added to modal inner products so the EVP
remains orthogonal under the chosen boundary law. Endpoint orientation
enters through Green's identity:
$$s_\ell=+1$$ at the surface and $$s_\ell=-1$$ at the bottom, so a
generated term with coefficient $$c$$ contributes
$$s_\ell c\,u_i(z_\ell)v_j(z_\ell)$$.

```matlab
boundaryConditions = IMBoundary.conditions(formulation="G", ...
    surface=IMBoundary.free(), bottom=IMBoundary.rigid());
```




## Topics
+ Create boundary conditions
  + [`active`](/internal-modes/classes-v2/core/imboundary/active.html) Create an active metadata-only boundary condition.
  + [`custom`](/internal-modes/classes-v2/core/imboundary/custom.html) Create a custom location-free operator boundary law.
  + [`dirichlet`](/internal-modes/classes-v2/core/imboundary/dirichlet.html) Create a location-free homogeneous Dirichlet boundary law.
  + [`free`](/internal-modes/classes-v2/core/imboundary/free.html) Create a location-free free boundary law.
  + [`innerProductTerm`](/internal-modes/classes-v2/core/imboundary/innerproductterm.html) Create a boundary inner-product trace-pair term.
  + [`linearF`](/internal-modes/classes-v2/core/imboundary/linearf.html) Create a location-free linear `F` boundary law.
  + [`linearG`](/internal-modes/classes-v2/core/imboundary/linearg.html) Create a location-free linear `G` boundary law.
  + [`neumann`](/internal-modes/classes-v2/core/imboundary/neumann.html) Create a location-free homogeneous Neumann boundary law.
  + [`noSlip`](/internal-modes/classes-v2/core/imboundary/noslip.html) Create a location-free no-slip boundary law.
  + [`partialDepthPE`](/internal-modes/classes-v2/core/imboundary/partialdepthpe.html) Create partial-depth potential-energy active boundary conditions.
  + [`rigid`](/internal-modes/classes-v2/core/imboundary/rigid.html) Create a location-free rigid boundary law.
  + [`trace`](/internal-modes/classes-v2/core/imboundary/trace.html) Create an endpoint trace descriptor.
+ Place boundary conditions
  + [`at`](/internal-modes/classes-v2/core/imboundary/at.html) Place a boundary condition at a physical endpoint.
  + [`conditions`](/internal-modes/classes-v2/core/imboundary/conditions.html) Place bottom and surface boundary conditions for one variable.
+ Inspect boundary conditions
  + [`boundaryModeDescriptors`](/internal-modes/classes-v2/core/imboundary/boundarymodedescriptors.html) Return declared endpoint boundary-mode metadata.
  + [`boundaryModeNumber`](/internal-modes/classes-v2/core/imboundary/boundarymodenumber.html) Physical mode number for an endpoint boundary branch.
  + [`expectedNegativeCount`](/internal-modes/classes-v2/core/imboundary/expectednegativecount.html) Return the negative-index contribution from this condition.
  + [`expectedZeroCount`](/internal-modes/classes-v2/core/imboundary/expectedzerocount.html) Return the zero-index contribution from this condition.
  + [`family`](/internal-modes/classes-v2/core/imboundary/family.html) Boundary family identifier.
  + [`hasInnerProductTerms`](/internal-modes/classes-v2/core/imboundary/hasinnerproductterms.html) Return true when this condition contributes inner-product terms.
  + [`hasKnownInnerProductTerms`](/internal-modes/classes-v2/core/imboundary/hasknowninnerproductterms.html) True when the compatible boundary inner-product terms are known.
  + [`indexRank`](/internal-modes/classes-v2/core/imboundary/indexrank.html) Rank of the boundary-mode eigenvalue contribution.
  + [`indexSign`](/internal-modes/classes-v2/core/imboundary/indexsign.html) Sign of the boundary-mode eigenvalue contribution.
  + [`innerProductTerms`](/internal-modes/classes-v2/core/imboundary/innerproductterms.html) Boundary terms implied by this condition for modal inner products.
  + [`leftOperator`](/internal-modes/classes-v2/core/imboundary/leftoperator.html) Left-side boundary functional.
  + [`location`](/internal-modes/classes-v2/core/imboundary/location.html) Physical endpoint location.
  + [`rightOperator`](/internal-modes/classes-v2/core/imboundary/rightoperator.html) Right-side eigenvalue boundary functional.
  + [`variable`](/internal-modes/classes-v2/core/imboundary/variable.html) constrained by this boundary condition.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`emptyInnerProductTerms`](/internal-modes/classes-v2/core/imboundary/emptyinnerproductterms.html) Return an empty boundary inner-product-term structure.


---