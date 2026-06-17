---
layout: default
title: IMAnalyticalInternalModesBasis
has_children: false
has_toc: false
mathjax: true
parent: Analytical bases
grand_parent: Class documentation V2
nav_order: 4
---

#  IMAnalyticalInternalModesBasis

Store exact internal-mode functions from an analytical solution family.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMAnalyticalInternalModesBasis</code></pre></div></div>

## Overview

`IMAnalyticalInternalModesBasis` has the same user-facing `F(z)`,
`G(z)`, normalization, and Gram-matrix idiom as numerical internal-mode
bases, but it evaluates closed-form functions supplied by an
`IMAnalyticalSolution`.

```matlab
solution = IMConstantStratificationSolution(N0=5.2e-3,zDomain=[-5000 0]);
basisSet = solution.internalModes(evp,nModes=4);
G = basisSet.G(linspace(-5000,0,128).');
```




## Topics
+ Evaluate analytical modes
  + [`F`](/internal-modes/classes-v2/analytical-bases/imanalyticalinternalmodesbasis/f.html) Evaluate exact `F` modes.
  + [`G`](/internal-modes/classes-v2/analytical-bases/imanalyticalinternalmodesbasis/g.html) Evaluate exact `G` modes.
  + [`IMAnalyticalInternalModesBasis`](/internal-modes/classes-v2/analytical-bases/imanalyticalinternalmodesbasis/imanalyticalinternalmodesbasis.html) Create an exact internal-mode basis.
  + [`addNormalization`](/internal-modes/classes-v2/analytical-bases/imanalyticalinternalmodesbasis/addnormalization.html) Add a named normalization rule.
  + [`normalization`](/internal-modes/classes-v2/analytical-bases/imanalyticalinternalmodesbasis/normalization.html) Active normalization rule name.
  + [`normalizationFactors`](/internal-modes/classes-v2/analytical-bases/imanalyticalinternalmodesbasis/normalizationfactors.html) Return scale factors for a normalization rule.
  + [`uz`](/internal-modes/classes-v2/analytical-bases/imanalyticalinternalmodesbasis/uz.html) Evaluate exact solved-variable derivatives.
+ Analyze Gram matrices
  + [`crossSpectrum`](/internal-modes/classes-v2/analytical-bases/imanalyticalinternalmodesbasis/crossspectrum.html) Compute a modal cross-spectrum.
  + [`gramMatrix`](/internal-modes/classes-v2/analytical-bases/imanalyticalinternalmodesbasis/grammatrix.html) Return a Gram matrix for exact `F` or `G` modes.
  + [`partialWindowModes`](/internal-modes/classes-v2/analytical-bases/imanalyticalinternalmodesbasis/partialwindowmodes.html) Diagonalize a partial-depth Gram matrix.
  + [`spectrum`](/internal-modes/classes-v2/analytical-bases/imanalyticalinternalmodesbasis/spectrum.html) Compute a modal spectrum.
+ Inspect analytical modes
  + [`N2`](/internal-modes/classes-v2/analytical-bases/imanalyticalinternalmodesbasis/n2.html) Buoyancy frequency squared function.
  + [`eigenvalues`](/internal-modes/classes-v2/analytical-bases/imanalyticalinternalmodesbasis/eigenvalues.html) Retained eigenvalues.
  + [`evp`](/internal-modes/classes-v2/analytical-bases/imanalyticalinternalmodesbasis/evp.html) Internal-mode EVP represented by these exact functions.
  + [`h`](/internal-modes/classes-v2/analytical-bases/imanalyticalinternalmodesbasis/h.html) Equivalent depths for retained modes.
  + [`metadata`](/internal-modes/classes-v2/analytical-bases/imanalyticalinternalmodesbasis/metadata.html) Additional creation metadata.
  + [`modeNumber`](/internal-modes/classes-v2/analytical-bases/imanalyticalinternalmodesbasis/modenumber.html) Retained mode labels.
  + [`modeSelectionDiagnostics`](/internal-modes/classes-v2/analytical-bases/imanalyticalinternalmodesbasis/modeselectiondiagnostics.html) Mode-selection diagnostics.
  + [`normalizationNames`](/internal-modes/classes-v2/analytical-bases/imanalyticalinternalmodesbasis/normalizationnames.html) Return installed normalization rule names.
  + [`solution`](/internal-modes/classes-v2/analytical-bases/imanalyticalinternalmodesbasis/solution.html) Analytical solution family that created this basis.
  + [`zDomain`](/internal-modes/classes-v2/analytical-bases/imanalyticalinternalmodesbasis/zdomain.html) Physical vertical domain.


---