---
layout: default
title: IMBoundaryWeight
has_children: false
has_toc: false
mathjax: true
parent: Core
grand_parent: Class documentation V2
nav_order: 5
---

#  IMBoundaryWeight

Describe one endpoint contribution to a modal inner product.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMBoundaryWeight</code></pre></div></div>

## Overview

`IMBoundaryWeight` is an atomic bilinear endpoint term. The
`innerProduct` property names the modal inner product receiving the
contribution, while the left and right endpoint factors identify which
modal variables and physical derivative orders are evaluated.

A resolved weight contributes
$$c\,X_i^{(p)}(z_\ell)Y_j^{(q)}(z_\ell)$$
at `location`, where `coefficient` is $$c$$, `leftVariable` is $$X$$,
`leftDerivativeOrder` is $$p$$, `rightVariable` is $$Y$$, and
`rightDerivativeOrder` is $$q$$. Location-free weights may be stored on
a boundary law and are placed when the law is resolved at the surface
or bottom endpoint.

```matlab
weight = IMBoundaryWeight(innerProduct="G", location="surface", ...
    coefficient=1, leftVariable="G", rightVariable="G");
```




## Topics
+ Create boundary weights
  + [`IMBoundaryWeight`](/internal-modes/classes-v2/core/imboundaryweight/imboundaryweight.html) Create an endpoint inner-product weight.
+ Inspect boundary weights
  + [`coefficient`](/internal-modes/classes-v2/core/imboundaryweight/coefficient.html) Endpoint coefficient.
  + [`innerProduct`](/internal-modes/classes-v2/core/imboundaryweight/innerproduct.html) Inner product receiving this endpoint contribution.
  + [`leftDerivativeOrder`](/internal-modes/classes-v2/core/imboundaryweight/leftderivativeorder.html) Physical derivative order of the left mode factor.
  + [`leftVariable`](/internal-modes/classes-v2/core/imboundaryweight/leftvariable.html) Variable evaluated for the left mode factor.
  + [`location`](/internal-modes/classes-v2/core/imboundaryweight/location.html) Physical endpoint location.
  + [`rightDerivativeOrder`](/internal-modes/classes-v2/core/imboundaryweight/rightderivativeorder.html) Physical derivative order of the right mode factor.
  + [`rightVariable`](/internal-modes/classes-v2/core/imboundaryweight/rightvariable.html) Variable evaluated for the right mode factor.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`at`](/internal-modes/classes-v2/core/imboundaryweight/at.html) Place location-free weights at a physical endpoint.
  + [`isScalarEndpointValue`](/internal-modes/classes-v2/core/imboundaryweight/isscalarendpointvalue.html) Return true for scalar same-variable endpoint-value weights.


---