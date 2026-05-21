---
layout: default
title: IMIndexPolicy
has_children: false
has_toc: false
mathjax: true
parent: Core
grand_parent: Class documentation V2
nav_order: 5
---

#  IMIndexPolicy

Specify the expected eigenvalue index of an EVP.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMIndexPolicy</code></pre></div></div>

## Overview

The index records negative, zero, and positive eigenvalue counts. Mode
numbers are separate physical labels: `-1` identifies a surface
boundary branch, `-2` identifies a bottom boundary branch, `0`
identifies a true null mode, and positive labels identify interior
baroclinic modes.

```matlab
boundaryConditions = IMBoundary.partialDepthPE(boundarySign="negative");
policy = IMIndexPolicy.fromBoundaryConditions(boundaryConditions);
```




## Topics
+ Create index policies
  + [`IMIndexPolicy`](/internal-modes/classes-v2/core/imindexpolicy/imindexpolicy.html) Create an index policy.
  + [`fixed`](/internal-modes/classes-v2/core/imindexpolicy/fixed.html) Create a policy with fixed expected negative and zero counts.
  + [`fromBoundaryConditions`](/internal-modes/classes-v2/core/imindexpolicy/fromboundaryconditions.html) Create an index policy from boundary-condition index metadata.
  + [`fromBoundarySigns`](/internal-modes/classes-v2/core/imindexpolicy/fromboundarysigns.html) Create an index policy from active-boundary signs.
  + [`none`](/internal-modes/classes-v2/core/imindexpolicy/none.html) Create a policy that records but does not validate index counts.
+ Validate index counts
  + [`boundaryModes`](/internal-modes/classes-v2/core/imindexpolicy/boundarymodes.html) Declared endpoint boundary-mode slots.
  + [`classify`](/internal-modes/classes-v2/core/imindexpolicy/classify.html) eigenvalues and validate the observed index.
  + [`expectedNegativeCount`](/internal-modes/classes-v2/core/imindexpolicy/expectednegativecount.html) Return the expected negative index count.
  + [`expectedNegativeCountValue`](/internal-modes/classes-v2/core/imindexpolicy/expectednegativecountvalue.html) Expected number of negative directions.
  + [`expectedPositiveCount`](/internal-modes/classes-v2/core/imindexpolicy/expectedpositivecount.html) Return the expected positive index count.
  + [`expectedZeroCount`](/internal-modes/classes-v2/core/imindexpolicy/expectedzerocount.html) Return the expected zero index count.
  + [`expectedZeroCountValue`](/internal-modes/classes-v2/core/imindexpolicy/expectedzerocountvalue.html) Expected number of zero directions.
  + [`indexTolerance`](/internal-modes/classes-v2/core/imindexpolicy/indextolerance.html) Eigenvalue tolerance for zero-index classification.
  + [`selectModes`](/internal-modes/classes-v2/core/imindexpolicy/selectmodes.html) Select and label retained modes according to the index policy.
  + [`validationMode`](/internal-modes/classes-v2/core/imindexpolicy/validationmode.html) Validation behavior.


---