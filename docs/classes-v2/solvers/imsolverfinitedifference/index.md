---
layout: default
title: IMSolverFiniteDifference
has_children: false
has_toc: false
mathjax: true
parent: Solvers
grand_parent: Class documentation V2
nav_order: 5
---

#  IMSolverFiniteDifference

Solve physical-coordinate EVPs on a supplied finite-difference grid.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMSolverFiniteDifference</code></pre></div></div>

## Overview

`IMSolverFiniteDifference` uses the user's physical `z`
grid as its native basis. It shares the `solveEVP` interface with
the spectral solvers but evaluates modes by interpolation.

```matlab
solver = IMSolverFiniteDifference(z=linspace(-1000,0,65).');
```




## Topics
+ Create solvers
  + [`IMSolverFiniteDifference`](/internal-modes/classes-v2/solvers/imsolverfinitedifference/imsolverfinitedifference.html) Create a finite-difference solver from a physical grid.
+ Inspect solvers
  + [`nEVP`](/internal-modes/classes-v2/solvers/imsolverfinitedifference/nevp.html) Number of native EVP values.
  + [`xNative`](/internal-modes/classes-v2/solvers/imsolverfinitedifference/xnative.html) Native coordinate alias for API consistency.
  + [`zDomain`](/internal-modes/classes-v2/solvers/imsolverfinitedifference/zdomain.html) Physical vertical domain.
  + [`zNative`](/internal-modes/classes-v2/solvers/imsolverfinitedifference/znative.html) Native finite-difference grid.
+ Evaluate native modes
  + [`xOfZ`](/internal-modes/classes-v2/solvers/imsolverfinitedifference/xofz.html) Return the native coordinate for a physical finite-difference grid.
  + [`zOfX`](/internal-modes/classes-v2/solvers/imsolverfinitedifference/zofx.html) Return the physical coordinate for a native finite-difference grid.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`T`](/internal-modes/classes-v2/solvers/imsolverfinitedifference/t.html) Native value matrix.
  + [`Tx`](/internal-modes/classes-v2/solvers/imsolverfinitedifference/tx.html) Native first-derivative matrix.
  + [`Txx`](/internal-modes/classes-v2/solvers/imsolverfinitedifference/txx.html) Native second-derivative matrix.


---