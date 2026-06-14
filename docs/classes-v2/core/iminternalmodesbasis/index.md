---
layout: default
title: IMInternalModesBasis
has_children: false
has_toc: false
mathjax: true
parent: Core
grand_parent: Class documentation V2
nav_order: 5
---

#  IMInternalModesBasis

Store solved internal-mode basis functions.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMInternalModesBasis < IMBasisSet</code></pre></div></div>

## Overview

`IMInternalModesBasis` evaluates the physical `F` and `G` variables
from the solved canonical scalar mode. If the EVP solves `G`, then
$$F_j=h_j G'_j.$$ If the EVP solves `F`, then
$$G_j=-gN^{-2}F'_j.$$




## Topics
+ Evaluate internal-mode bases
  + [`F`](/internal-modes/classes-v2/core/iminternalmodesbasis/f.html) Evaluate `F` modes.
  + [`G`](/internal-modes/classes-v2/core/iminternalmodesbasis/g.html) Evaluate `G` modes.
  + [`IMInternalModesBasis`](/internal-modes/classes-v2/core/iminternalmodesbasis/iminternalmodesbasis.html) Create an internal-mode basis set.


---