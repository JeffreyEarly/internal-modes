---
layout: default
title: InternalModesAdaptiveSpectral
has_children: false
has_toc: false
mathjax: true
parent: Numerical solvers
grand_parent: Class documentation
nav_order: 3
---

#  InternalModesAdaptiveSpectral

Solve the vertical EVP on an adaptive WKB grid with coupled Chebyshev blocks.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef InternalModesAdaptiveSpectral < InternalModesWKBSpectral</code></pre></div></div>

## Overview

`InternalModesAdaptiveSpectral` extends
`InternalModesWKBSpectral` using the adaptive multi-region strategy
described in Section 4.4 of Early, Lelong, and Smith (2020). The
solver keeps the WKB stretched coordinate but partitions it into
oscillatory and evanescent regions, then couples separate Chebyshev
blocks across the turning points.

This is most useful at fixed frequency, where the turning points of
$$N^2(z) - \omega^2$$ can leave large regions that are exponentially
decaying rather than oscillatory.

```matlab
im = InternalModesAdaptiveSpectral(rho=rho, zIn=zIn, zOut=zOut, latitude=latitude, nEVP=257);
[F, G, h, k] = im.ModesAtFrequency(2*pi*1e-3);
```




## Topics
+ Create and initialize modes
  + [`InternalModesAdaptiveSpectral`](/internal-modes/classes/numerical-solvers/internalmodesadaptivespectral/internalmodesadaptivespectral.html) Initialize the adaptive WKB spectral solver.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Developer topics
  + [`Diff1_xChebFunction`](/internal-modes/classes/numerical-solvers/internalmodesadaptivespectral/diff1_xchebfunction.html) Differentiate a vector in the coupled Chebyshev basis.
  + [`DistributePointsInRegionsWithMinimum`](/internal-modes/classes/numerical-solvers/internalmodesadaptivespectral/distributepointsinregionswithminimum.html) Distribute collocation points across adaptive regions.
  + [`FindWKBSolutionBoundaries`](/internal-modes/classes/numerical-solvers/internalmodesadaptivespectral/findwkbsolutionboundaries.html) Estimate adaptive-region boundaries from turning points and WKB decay.
  + [`Lxi`](/internal-modes/classes/numerical-solvers/internalmodesadaptivespectral/lxi.html) Length of each coupled subdomain in the adaptive coordinate.
  + [`N_zCheb`](/internal-modes/classes/numerical-solvers/internalmodesadaptivespectral/n_zcheb.html) Chebyshev coefficients of $$N(z)$$ on the reference grid.
  + [`RemoveRegionAtIndex`](/internal-modes/classes/numerical-solvers/internalmodesadaptivespectral/removeregionatindex.html) Merge neighboring adaptive regions by removing one boundary.
  + [`T_xCheb_xLobatto`](/internal-modes/classes/numerical-solvers/internalmodesadaptivespectral/t_xcheb_xlobatto.html) Transform coupled Chebyshev coefficients onto the adaptive Lobatto grid.
  + [`T_xCheb_zOutFunction`](/internal-modes/classes/numerical-solvers/internalmodesadaptivespectral/t_xcheb_zoutfunction.html) Transform coupled Chebyshev coefficients onto the public output grid.
  + [`T_xCheb_zOut_Transforms`](/internal-modes/classes/numerical-solvers/internalmodesadaptivespectral/t_xcheb_zout_transforms.html) Per-region transforms from Chebyshev coefficients to `zOut`.
  + [`T_xCheb_zOut_fromIndices`](/internal-modes/classes/numerical-solvers/internalmodesadaptivespectral/t_xcheb_zout_fromindices.html) Source indices in coefficient space for each output transform.
  + [`T_xCheb_zOut_toIndices`](/internal-modes/classes/numerical-solvers/internalmodesadaptivespectral/t_xcheb_zout_toindices.html) Target indices in `zOut` for each output transform.
  + [`T_xLobatto_xCheb`](/internal-modes/classes/numerical-solvers/internalmodesadaptivespectral/t_xlobatto_xcheb.html) Transform adaptive-grid values into the coupled Chebyshev basis.
  + [`WKBDecaySolution`](/internal-modes/classes/numerical-solvers/internalmodesadaptivespectral/wkbdecaysolution.html) Evaluate the Airy-based WKB decay envelope used for pruning regions.
  + [`eqIndices`](/internal-modes/classes/numerical-solvers/internalmodesadaptivespectral/eqindices.html) Row indices for each coupled subproblem.
  + [`nEquations`](/internal-modes/classes/numerical-solvers/internalmodesadaptivespectral/nequations.html) Number of coupled Chebyshev subproblems.
  + [`polyIndices`](/internal-modes/classes/numerical-solvers/internalmodesadaptivespectral/polyindices.html) Column indices for each coupled Chebyshev block.
  + [`x_zLobatto`](/internal-modes/classes/numerical-solvers/internalmodesadaptivespectral/x_zlobatto.html) Adaptive stretched-coordinate grid sampled on the public Lobatto grid.
  + [`xiBoundaries`](/internal-modes/classes/numerical-solvers/internalmodesadaptivespectral/xiboundaries.html) WKB-coordinate values of region boundaries and turning points.
  + [`xiIndices`](/internal-modes/classes/numerical-solvers/internalmodesadaptivespectral/xiindices.html) Grid-point indices into the unique adaptive Lobatto grid.
  + [`zBoundaries`](/internal-modes/classes/numerical-solvers/internalmodesadaptivespectral/zboundaries.html) Physical depths of region boundaries and turning points.


---