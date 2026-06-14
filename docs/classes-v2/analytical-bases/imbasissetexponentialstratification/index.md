---
layout: default
title: IMBasisSetExponentialStratification
has_children: false
has_toc: false
mathjax: true
parent: Analytical bases
grand_parent: Class documentation V2
nav_order: 2
---

#  IMBasisSetExponentialStratification

Evaluate exact basis sets for exponential stratification.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMBasisSetExponentialStratification < IMInternalModesBasis</code></pre></div></div>

## Overview

`IMBasisSetExponentialStratification` stores exact rigid-bottom
`G`-formulation basis sets for rigid or free surfaces with
$$N^2(z)=N_0^2 e^{2z/b}$$ on domains with surface at $$z=0$$. Mode
roots are found with a local scanner and `fzero`, without requiring
an external spectral-function root finder at runtime.

```matlab
zDomain = [-5000 0];
N2 = @(z) (5.2e-3)^2*exp(2*z/1300);
evp = IMInternalModes.hydrostaticGModes(N2=N2, zDomain=zDomain);
basisSet = IMBasisSetExponentialStratification(evp=evp, N0=5.2e-3, b=1300, zDomain=[-5000 0]);
G = basisSet.G(linspace(-5000,0,128).');
```




## Topics
+ Create basis sets
  + [`IMBasisSetExponentialStratification`](/internal-modes/classes-v2/analytical-bases/imbasissetexponentialstratification/imbasissetexponentialstratification.html) Create an exact exponential-stratification basis set.
+ Inspect basis sets
  + [`N0`](/internal-modes/classes-v2/analytical-bases/imbasissetexponentialstratification/n0.html) Surface buoyancy frequency $$N_0$$ in radians per second.
  + [`b`](/internal-modes/classes-v2/analytical-bases/imbasissetexponentialstratification/b.html) Exponential e-folding depth $$b$$ in meters.
  + [`frequencies`](/internal-modes/classes-v2/analytical-bases/imbasissetexponentialstratification/frequencies.html) Modal frequencies in radians per second.
  + [`modeKinds`](/internal-modes/classes-v2/analytical-bases/imbasissetexponentialstratification/modekinds.html) Internal analytical branch kind for each retained mode.
  + [`phaseSpeeds`](/internal-modes/classes-v2/analytical-bases/imbasissetexponentialstratification/phasespeeds.html) Phase speeds $$c_j=\sqrt{g h_j}$$ for each retained mode.
  + [`roots`](/internal-modes/classes-v2/analytical-bases/imbasissetexponentialstratification/roots.html) Dimensionless Bessel roots used to construct the retained modes.
  + [`signFactors`](/internal-modes/classes-v2/analytical-bases/imbasissetexponentialstratification/signfactors.html) Sign applied to each raw mode so that $$F_j(0)>0$$.


---