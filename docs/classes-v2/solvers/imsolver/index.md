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

Define the shared protocol for canonical EVP solvers.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef (Abstract) IMSolver</code></pre></div></div>

## Overview

Concrete solvers own the grid, coordinate mapping, derivative
matrices, integration rule, and interpolation of native modes. The
base class owns the common generalized-eigenvalue workflow.




## Topics
+ Solve EVPs
  + [`nativeQuadratureRule`](/internal-modes/classes-v2/solvers/imsolver/nativequadraturerule.html) Return the solver's native physical quadrature rule.
  + [`solveEVP`](/internal-modes/classes-v2/solvers/imsolver/solveevp.html) Solve an EVP and return a basis set.
+ Solve geostrophic zero-APV modes
  + [`solveGeostrophicZeroAPVModes`](/internal-modes/classes-v2/solvers/imsolver/solvegeostrophiczeroapvmodes.html) Solve canonical geostrophic zero-APV boundary modes.
+ Other
  + [`IMSolver`](/internal-modes/classes-v2/solvers/imsolver/imsolver.html)
  + [`N2`](/internal-modes/classes-v2/solvers/imsolver/n2.html)
  + [`boundaryIndex`](/internal-modes/classes-v2/solvers/imsolver/boundaryindex.html)
  + [`configuredForEVP`](/internal-modes/classes-v2/solvers/imsolver/configuredforevp.html)
  + [`context`](/internal-modes/classes-v2/solvers/imsolver/context.html)
  + [`differentiateGridValues`](/internal-modes/classes-v2/solvers/imsolver/differentiategridvalues.html)
  + [`evaluateNativeModes`](/internal-modes/classes-v2/solvers/imsolver/evaluatenativemodes.html)
  + [`evaluatePhysicalDerivative`](/internal-modes/classes-v2/solvers/imsolver/evaluatephysicalderivative.html)
  + [`innerProductGrid`](/internal-modes/classes-v2/solvers/imsolver/innerproductgrid.html)
  + [`integrateGridValuesFromSurface`](/internal-modes/classes-v2/solvers/imsolver/integrategridvaluesfromsurface.html)
  + [`integrateInnerProduct`](/internal-modes/classes-v2/solvers/imsolver/integrateinnerproduct.html)
  + [`physicalDerivativeMatrix`](/internal-modes/classes-v2/solvers/imsolver/physicalderivativematrix.html)


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Solve geostrophic zero-APV modes
  + [`configuredForGeostrophicZeroAPVModes`](/internal-modes/classes-v2/solvers/imsolver/configuredforgeostrophiczeroapvmodes.html) Return a solver configured for geostrophic zero-APV modes.
+ Developer topics
  + [`rootsOfNativeMode`](/internal-modes/classes-v2/solvers/imsolver/rootsofnativemode.html) Return physical roots of one native mode.


---