---
layout: default
title: IMProductAssessment
has_children: false
has_toc: false
mathjax: true
parent: Discrete transforms
grand_parent: Class documentation V2
nav_order: 7
---

#  IMProductAssessment

Report sampled product errors separately from acceptance policy.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMProductAssessment</code></pre></div></div>

## Overview

Each measurement is a maximum over individually examined products.
Neither a passing sparse plan nor a bounded all-products control proves
arbitrary-superposition accuracy or complete physical qualification.

```matlab
result = plan.assess(grids);
decision = result.applyPolicy(quadraticTolerance=0.1,referenceTolerance=1e-4);
```




## Topics
+ Inspect product evidence
  + [`IMProductAssessment`](/internal-modes/classes-v2/discrete-transforms/improductassessment/improductassessment.html) Aggregate complete executed evidence without applying tolerances.
  + [`assessmentIdentity`](/internal-modes/classes-v2/discrete-transforms/improductassessment/assessmentidentity.html) Exact grids and prepared numerical output/reference identities.
  + [`compareCoverage`](/internal-modes/classes-v2/discrete-transforms/improductassessment/comparecoverage.html) Detect failures missed by this plan using a separate bounded control.
  + [`costs`](/internal-modes/classes-v2/discrete-transforms/improductassessment/costs.html) Planning, preparation, evaluation, assessment, and array memory costs.
  + [`coverage`](/internal-modes/classes-v2/discrete-transforms/improductassessment/coverage.html) Examined and omitted inventory families and interactions.
  + [`evidence`](/internal-modes/classes-v2/discrete-transforms/improductassessment/evidence.html) Individual products, signs, labels, zeros, and reference evidence.
  + [`measurements`](/internal-modes/classes-v2/discrete-transforms/improductassessment/measurements.html) One error/status/limiting-identity row per requested count.
  + [`plan`](/internal-modes/classes-v2/discrete-transforms/improductassessment/plan.html) Immutable reservation and original scientific inventory.
+ Apply product policies
  + [`applyPolicy`](/internal-modes/classes-v2/discrete-transforms/improductassessment/applypolicy.html) Apply explicit product/reference tolerances without reducing counts.


---