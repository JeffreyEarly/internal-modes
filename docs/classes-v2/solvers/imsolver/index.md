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
  + [`IMSolver`](/internal-modes/classes-v2/solvers/imsolver/imsolver.html) Define the shared protocol for canonical EVP solvers.
  + [`solveEVP`](/internal-modes/classes-v2/solvers/imsolver/solveevp.html) Solve an EVP and return a basis set.
+ Solve surface-geostrophic modes
  + [`solveSurfaceGeostrophicModes`](/internal-modes/classes-v2/solvers/imsolver/solvesurfacegeostrophicmodes.html) Solve projected surface-geostrophic boundary modes.
+ Other
  + [`N2`](/internal-modes/classes-v2/solvers/imsolver/n2.html)
  + [`boundaryIndex`](/internal-modes/classes-v2/solvers/imsolver/boundaryindex.html)
  + [`configuredForEVP`](/internal-modes/classes-v2/solvers/imsolver/configuredforevp.html)
  + [`context`](/internal-modes/classes-v2/solvers/imsolver/context.html)
  + [`differentiateGridValues`](/internal-modes/classes-v2/solvers/imsolver/differentiategridvalues.html)
  + [`evaluateNativeModes`](/internal-modes/classes-v2/solvers/imsolver/evaluatenativemodes.html)
  + [`evaluatePhysicalDerivative`](/internal-modes/classes-v2/solvers/imsolver/evaluatephysicalderivative.html)
  + [`innerProductGrid`](/internal-modes/classes-v2/solvers/imsolver/innerproductgrid.html)
  + [`integrateInnerProduct`](/internal-modes/classes-v2/solvers/imsolver/integrateinnerproduct.html)
  + [`physicalDerivativeMatrix`](/internal-modes/classes-v2/solvers/imsolver/physicalderivativematrix.html)


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Solve surface-geostrophic modes
  + [`configuredForSurfaceGeostrophicModes`](/internal-modes/classes-v2/solvers/imsolver/configuredforsurfacegeostrophicmodes.html) Return a solver configured for surface-geostrophic modes.


---