---
layout: default
title: IMSolver
has_children: false
has_toc: false
mathjax: true
parent: Solvers
grand_parent: Class documentation V2
nav_order: 1
---

#  IMSolver

Define the shared protocol for internal-mode solvers.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef (Abstract) IMSolver</code></pre></div></div>

## Overview

`IMSolver` owns the solver-independent generalized EVP
workflow. Concrete subclasses provide the native grid, physical
derivative matrices, boundary rows, and native-mode evaluation.
Solvers own the numerical medium and discretization. EVPs own the
physical constants and combine them with solver context during
assembly.




## Topics
+ Solve EVPs
  + [`IMSolver`](/internal-modes/classes-v2/solvers/imsolver/imsolver.html) Define the shared protocol for internal-mode solvers.
  + [`solveEVP`](/internal-modes/classes-v2/solvers/imsolver/solveevp.html) Solve an EVP and return a native-basis solution set.
+ Other
  + [`N2`](/internal-modes/classes-v2/solvers/imsolver/n2.html)
  + [`boundaryIndex`](/internal-modes/classes-v2/solvers/imsolver/boundaryindex.html)
  + [`context`](/internal-modes/classes-v2/solvers/imsolver/context.html)
  + [`dzLogN2`](/internal-modes/classes-v2/solvers/imsolver/dzlogn2.html)
  + [`evaluateNativeModes`](/internal-modes/classes-v2/solvers/imsolver/evaluatenativemodes.html)
  + [`evaluatePhysicalDerivative`](/internal-modes/classes-v2/solvers/imsolver/evaluatephysicalderivative.html)
  + [`innerProductGrid`](/internal-modes/classes-v2/solvers/imsolver/innerproductgrid.html)
  + [`integrateInnerProduct`](/internal-modes/classes-v2/solvers/imsolver/integrateinnerproduct.html)
  + [`physicalDerivativeMatrix`](/internal-modes/classes-v2/solvers/imsolver/physicalderivativematrix.html)


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`applyEndpointLaw`](/internal-modes/classes-v2/solvers/imsolver/applyendpointlaw.html) Apply a resolved endpoint law to a matrix pair.


---