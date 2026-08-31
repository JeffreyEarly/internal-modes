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

Describe common state for a stratification family with exact solutions.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef (Abstract) IMAnalyticalSolution</code></pre></div></div>

## Overview

`IMAnalyticalSolution` owns only physical state shared by analytical
stratification families. Concrete classes advertise exact solution
operations by implementing those construction methods directly.




## Topics
+ Create analytical solutions
  + [`IMAnalyticalSolution`](/internal-modes/classes-v2/analytical-bases/imanalyticalsolution/imanalyticalsolution.html) Initialize common analytical-stratification state.
+ Inspect analytical solutions
  + [`N2`](/internal-modes/classes-v2/analytical-bases/imanalyticalsolution/n2.html) Evaluate buoyancy frequency squared $$N^2(z)$$.
  + [`f0`](/internal-modes/classes-v2/analytical-bases/imanalyticalsolution/f0.html) Coriolis parameter in radians per second.
  + [`g`](/internal-modes/classes-v2/analytical-bases/imanalyticalsolution/g.html) Gravitational acceleration in meters per second squared.
  + [`summarize`](/internal-modes/classes-v2/analytical-bases/imanalyticalsolution/summarize.html) Print common analytical-stratification state.
  + [`zDomain`](/internal-modes/classes-v2/analytical-bases/imanalyticalsolution/zdomain.html) Physical vertical domain.


---