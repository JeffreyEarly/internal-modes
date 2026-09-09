---
layout: default
title: IMProductInventory
has_children: false
has_toc: false
mathjax: true
parent: Discrete transforms
grand_parent: Class documentation V2
nav_order: 5
---

#  IMProductInventory

Declare products, continuous factors, and caller-owned output recipes.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMProductInventory</code></pre></div></div>

## Overview

Construction inspects metadata only. Physical channels, valid
interactions, frequency signs, and endpoint recipes belong to the
caller. Factors and output recipes are evaluated only after a plan
reserves every selected individual product, including exact zeros.

```matlab
inventory = IMProductInventory(factors,products,outputs);
plan = inventory.fixedPlan(prefixCounts=1:8,productBudget=10000);
```




## Topics
+ Declare inventories
  + [`IMProductInventory`](/internal-modes/classes-v2/discrete-transforms/improductinventory/improductinventory.html) Validate an inventory without invoking any evaluator.
  + [`collectionFactor`](/internal-modes/classes-v2/discrete-transforms/improductinventory/collectionfactor.html) Bind one collection page, variable, derivative, and coefficient function.
  + [`factors`](/internal-modes/classes-v2/discrete-transforms/improductinventory/factors.html) Cell row of factor metadata and explicit evaluation callbacks.
  + [`outputs`](/internal-modes/classes-v2/discrete-transforms/improductinventory/outputs.html) Cell row of metadata and deferred output preparation callbacks.
  + [`products`](/internal-modes/classes-v2/discrete-transforms/improductinventory/products.html) Ordered table of valid interaction/channel/factor/output identities.
+ Plan assessments
  + [`allProductsPlan`](/internal-modes/classes-v2/discrete-transforms/improductinventory/allproductsplan.html) Reserve all supplied products for bounded independent validation.
  + [`fixedPlan`](/internal-modes/classes-v2/discrete-transforms/improductinventory/fixedplan.html) Reserve a deterministic sparse plan before numerical preparation.


---