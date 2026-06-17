---
layout: default
title: IMAnalyticalSolution
has_children: false
has_toc: false
mathjax: true
parent: Analytical bases
grand_parent: Class documentation V2
nav_order: 1
---

#  IMAnalyticalSolution

Describe a stratification family with closed-form mode solutions.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMAnalyticalSolution</code></pre></div></div>

## Overview

`IMAnalyticalSolution` is the base class for exact stratification
families. A concrete solution object owns the continuous stratification
parameters and advertises which internal-mode and SQG branches are
available.

```matlab
solution = IMConstantStratificationSolution(N0=5.2e-3,zDomain=[-5000 0]);
availability = solution.internalModeAvailability(evp);
```




## Topics
+ Create analytical solutions
  + [`IMAnalyticalSolution`](/internal-modes/classes-v2/analytical-bases/imanalyticalsolution/imanalyticalsolution.html) Create an analytical solution family.
+ Inspect analytical solutions
  + [`N2`](/internal-modes/classes-v2/analytical-bases/imanalyticalsolution/n2.html) Evaluate the stratification profile.
  + [`f0`](/internal-modes/classes-v2/analytical-bases/imanalyticalsolution/f0.html) Coriolis parameter in radians per second.
  + [`g`](/internal-modes/classes-v2/analytical-bases/imanalyticalsolution/g.html) Gravitational acceleration in meters per second squared.
  + [`summarize`](/internal-modes/classes-v2/analytical-bases/imanalyticalsolution/summarize.html) Print a readable summary of the analytical solution.
  + [`zDomain`](/internal-modes/classes-v2/analytical-bases/imanalyticalsolution/zdomain.html) Physical vertical domain.
+ Compute internal modes
  + [`internalModeAvailability`](/internal-modes/classes-v2/analytical-bases/imanalyticalsolution/internalmodeavailability.html) Report whether exact internal modes are available.
  + [`internalModes`](/internal-modes/classes-v2/analytical-bases/imanalyticalsolution/internalmodes.html) Create an exact internal-mode basis when available.
+ Compute SQG modes
  + [`sqgAvailability`](/internal-modes/classes-v2/analytical-bases/imanalyticalsolution/sqgavailability.html) Report whether exact SQG boundary modes are available.
  + [`sqgModesAtWavenumber`](/internal-modes/classes-v2/analytical-bases/imanalyticalsolution/sqgmodesatwavenumber.html) Create exact SQG boundary modes at fixed wavenumber.


---