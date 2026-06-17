---
layout: default
title: IMAnalyticalSQGBasis
has_children: false
has_toc: false
mathjax: true
parent: Analytical bases
grand_parent: Class documentation V2
nav_order: 5
---

#  IMAnalyticalSQGBasis

Store exact SQG boundary modes from an analytical solution family.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMAnalyticalSQGBasis</code></pre></div></div>

## Overview

`IMAnalyticalSQGBasis` evaluates the boundary-trapped SQG streamfunction
modes available for a closed-form stratification family.

```matlab
solution = IMConstantStratificationSolution(N0=5.2e-3,zDomain=[-5000 0]);
sqg = solution.sqgModesAtWavenumber(1e-4,boundary="surface");
psi = sqg.psi(z);
```




## Topics
+ Evaluate SQG modes
  + [`IMAnalyticalSQGBasis`](/internal-modes/classes-v2/analytical-bases/imanalyticalsqgbasis/imanalyticalsqgbasis.html) Create an exact SQG basis.
  + [`psi`](/internal-modes/classes-v2/analytical-bases/imanalyticalsqgbasis/psi.html) Evaluate exact SQG streamfunction modes.
+ Inspect SQG modes
  + [`N2`](/internal-modes/classes-v2/analytical-bases/imanalyticalsqgbasis/n2.html) Buoyancy frequency squared function.
  + [`boundary`](/internal-modes/classes-v2/analytical-bases/imanalyticalsqgbasis/boundary.html) Active SQG boundary, `"surface"` or `"bottom"`.
  + [`k`](/internal-modes/classes-v2/analytical-bases/imanalyticalsqgbasis/k.html) Horizontal wavenumbers.
  + [`metadata`](/internal-modes/classes-v2/analytical-bases/imanalyticalsqgbasis/metadata.html) Additional creation metadata.
  + [`solution`](/internal-modes/classes-v2/analytical-bases/imanalyticalsqgbasis/solution.html) Analytical solution family that created this SQG basis.
  + [`summarize`](/internal-modes/classes-v2/analytical-bases/imanalyticalsqgbasis/summarize.html) Print a readable SQG basis summary.
  + [`zDomain`](/internal-modes/classes-v2/analytical-bases/imanalyticalsqgbasis/zdomain.html) Physical vertical domain.


---