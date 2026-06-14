---
layout: default
title: IMInternalModes
has_children: false
has_toc: false
mathjax: true
parent: Core
grand_parent: Class documentation V2
nav_order: 4
---

#  IMInternalModes

Describe canonical EVPs with internal-mode interpretation.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMInternalModes < IMEigenvalueProblem</code></pre></div></div>

## Overview

`IMInternalModes` translates standard `F` and `G` internal-mode
problems into the canonical scalar EVP. The solved scalar `u` is
either `F` or `G`; the other variable is recovered diagnostically by
the resulting `IMInternalModesBasis`.




## Topics
+ Create internal-mode EVPs
  + [`IMInternalModes`](/internal-modes/classes-v2/core/iminternalmodes/iminternalmodes.html) Create an internal-mode canonical EVP.
  + [`hydrostaticFModes`](/internal-modes/classes-v2/core/iminternalmodes/hydrostaticfmodes.html) Create the hydrostatic `F` internal-mode EVP.
  + [`hydrostaticGModes`](/internal-modes/classes-v2/core/iminternalmodes/hydrostaticgmodes.html) Create the hydrostatic `G` internal-mode EVP.
  + [`waveModesAtFrequency`](/internal-modes/classes-v2/core/iminternalmodes/wavemodesatfrequency.html) Create the fixed-frequency wave-mode EVP.
  + [`waveModesAtWavenumber`](/internal-modes/classes-v2/core/iminternalmodes/wavemodesatwavenumber.html) Create the fixed-wavenumber wave-mode EVP.
+ Inspect internal-mode EVPs
  + [`f0`](/internal-modes/classes-v2/core/iminternalmodes/f0.html) Coriolis parameter.
  + [`formulation`](/internal-modes/classes-v2/core/iminternalmodes/formulation.html) Solved physical variable, `"F"` or `"G"`.
  + [`g`](/internal-modes/classes-v2/core/iminternalmodes/g.html) Gravitational acceleration.


---