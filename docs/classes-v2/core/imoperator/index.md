---
layout: default
title: IMOperator
has_children: false
has_toc: false
mathjax: true
parent: Core
grand_parent: Class documentation V2
nav_order: 3
---

#  IMOperator

Represent a physical-coordinate linear differential operator.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMOperator</code></pre></div></div>

## Overview

`IMOperator` stores structured terms of the form
$$a_p(z)\partial_z^p$$. Coordinate-aware solvers pull these terms
back to their native coordinate when assembling EVP matrices.
Coefficient handles receive the EVP-built framework context `ctx`,
including `ctx.N2(z)`, `ctx.dzLogN2(z)`, `ctx.g`, `ctx.f0`,
`ctx.zDomain`, and `ctx.coordinateKind`.

```matlab
op = IMOperator().plus(coefficient=@(z,ctx) ctx.N2(z), derivativeOrder=0);
```




## Topics
+ Create operators
  + [`IMOperator`](/internal-modes/classes-v2/core/imoperator/imoperator.html) Create an empty physical-coordinate operator.
  + [`form`](/internal-modes/classes-v2/core/imoperator/form.html) Operator representation.
  + [`plus`](/internal-modes/classes-v2/core/imoperator/plus.html) Add a structured term to the operator.
+ Inspect operators
  + [`terms`](/internal-modes/classes-v2/core/imoperator/terms.html) Structured operator terms.
+ Assemble operators
  + [`boundaryRow`](/internal-modes/classes-v2/core/imoperator/boundaryrow.html) Assemble a boundary functional row.
  + [`matrix`](/internal-modes/classes-v2/core/imoperator/matrix.html) Assemble the operator on a solver's native basis.
+ Other
  + [`resolveContext`](/internal-modes/classes-v2/core/imoperator/resolvecontext.html)


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Create operators
  + [`strong`](/internal-modes/classes-v2/core/imoperator/strong.html) Create a strong-form physical-coordinate operator.
+ Assemble operators
  + [`evaluate`](/internal-modes/classes-v2/core/imoperator/evaluate.html) an operator applied to native mode columns.
  + [`evaluateCoefficient`](/internal-modes/classes-v2/core/imoperator/evaluatecoefficient.html) Evaluate an operator coefficient on a grid.


---