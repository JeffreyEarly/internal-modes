---
layout: default
title: IMProjectionRecipe
has_children: false
has_toc: false
mathjax: true
parent: Discrete transforms
grand_parent: Class documentation V2
nav_order: 2
---

#  IMProjectionRecipe

Bind continuous basis mathematics to fixed-grid numerical projection.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMProjectionRecipe</code></pre></div></div>

## Overview

The signed sample metric defines projection; the positive majorant
defines error magnitudes. This value object snapshots its basis and
normalization. It does not choose a mode count or certify solve accuracy.

```matlab
recipe = basis.projectionRecipe(variable="F");
projection = recipe.projection(z,weights);
```




## Topics
+ Create projection recipes
  + [`IMProjectionRecipe`](/internal-modes/classes-v2/discrete-transforms/improjectionrecipe/improjectionrecipe.html) Store a basis-owned recipe without changing its scientific state.
+ Inspect projection recipes
  + [`available`](/internal-modes/classes-v2/discrete-transforms/improjectionrecipe/available.html) Whether the continuous projection pairing is available.
  + [`columnLabels`](/internal-modes/classes-v2/discrete-transforms/improjectionrecipe/columnlabels.html) Scientific labels in basis-column order.
  + [`majorantGramMatrix`](/internal-modes/classes-v2/discrete-transforms/improjectionrecipe/majorantgrammatrix.html) Continuous positive-majorant Gram matrix.
  + [`normalization`](/internal-modes/classes-v2/discrete-transforms/improjectionrecipe/normalization.html) Snapshotted normalization convention.
  + [`provenance`](/internal-modes/classes-v2/discrete-transforms/improjectionrecipe/provenance.html) Representation and reference-integration provenance.
  + [`reason`](/internal-modes/classes-v2/discrete-transforms/improjectionrecipe/reason.html) Explanation when the continuous pairing is unavailable.
  + [`supportsLeakage`](/internal-modes/classes-v2/discrete-transforms/improjectionrecipe/supportsleakage.html) Whether the family supports leakage measurements.
  + [`supportsQuadratic`](/internal-modes/classes-v2/discrete-transforms/improjectionrecipe/supportsquadratic.html) Whether the family supports quadratic measurements.
  + [`targetGramMatrix`](/internal-modes/classes-v2/discrete-transforms/improjectionrecipe/targetgrammatrix.html) Continuous signed target Gram matrix.
  + [`variable`](/internal-modes/classes-v2/discrete-transforms/improjectionrecipe/variable.html) Scalar variable evaluated by this recipe.
+ Evaluate projection recipes
  + [`evaluate`](/internal-modes/classes-v2/discrete-transforms/improjectionrecipe/evaluate.html) normalized continuous columns at physical points.
  + [`metric`](/internal-modes/classes-v2/discrete-transforms/improjectionrecipe/metric.html) Construct the signed or majorant pairing on fixed samples.
  + [`projection`](/internal-modes/classes-v2/discrete-transforms/improjectionrecipe/projection.html) Construct a projection for explicit columns on an unchanged grid.


---