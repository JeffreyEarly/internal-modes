---
layout: default
title: IMBasisSetConstantStratification
has_children: false
has_toc: false
mathjax: true
parent: Analytical bases
grand_parent: Class documentation V2
nav_order: 1
---

#  IMBasisSetConstantStratification

Evaluate exact basis sets for constant stratification.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMBasisSetConstantStratification < IMBasisSet</code></pre></div></div>

## Overview

`IMBasisSetConstantStratification` stores exact
constant-stratification basis sets for $$N^2(z)=N_0^2$$. The class
implements the same `F`, `G`, and normalization contract as numerical
basis sets, without storing a solver reference.

```matlab
evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4);
basisSet = IMBasisSetConstantStratification(evp=evp, N0=5.2e-3, zDomain=[-5000 0]);
G = basisSet.G(linspace(-5000,0,128).');
```




## Topics
+ Create basis sets
  + [`IMBasisSetConstantStratification`](/internal-modes/classes-v2/analytical-bases/imbasissetconstantstratification/imbasissetconstantstratification.html) Create an exact constant-stratification basis set.
+ Inspect basis sets
  + [`N0`](/internal-modes/classes-v2/analytical-bases/imbasissetconstantstratification/n0.html) Constant buoyancy frequency $$N_0$$ in radians per second.
  + [`isBoundaryMode`](/internal-modes/classes-v2/analytical-bases/imbasissetconstantstratification/isboundarymode.html) True for a free-surface boundary branch.
  + [`solutionTypes`](/internal-modes/classes-v2/analytical-bases/imbasissetconstantstratification/solutiontypes.html) Analytical branch type for each retained mode.
  + [`verticalWavenumbers`](/internal-modes/classes-v2/analytical-bases/imbasissetconstantstratification/verticalwavenumbers.html) Vertical wavenumbers for each retained mode.


---