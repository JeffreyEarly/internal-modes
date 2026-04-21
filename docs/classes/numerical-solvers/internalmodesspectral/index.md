---
layout: default
title: InternalModesSpectral
has_children: false
has_toc: false
mathjax: true
parent: Numerical solvers
grand_parent: Class documentation
nav_order: 1
---

#  InternalModesSpectral

Solve the vertical internal-mode EVP with Chebyshev collocation.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef InternalModesSpectral < InternalModesBase</code></pre></div></div>

## Overview

`InternalModesSpectral` is the core numerical solver described in
Sections 4 and 5 of Early, Lelong, and Smith (2020). It represents
the vertical structure on a Gauss-Lobatto grid and solves the
manuscript eigenvalue problems

$$
\partial_{zz} G_j - K^2 G_j = \frac{f_0^2 - N^2}{g h_j} G_j
$$

for fixed `K`, or

$$
\partial_{zz} G_j = \frac{\omega^2 - N^2}{g h_j} G_j
$$

for fixed `\omega`, together with

$$
F_j = h_j \partial_z G_j.
$$

The solver stores background fields on the user-facing output grid
`z`, but internally works with a spectral grid and exact integrals
of Chebyshev polynomials for mode normalization. Subclasses reuse
this machinery in stretched coordinates while preserving the same
public mode API.

```matlab
im = InternalModesSpectral(rho=rho, zIn=zIn, zOut=zOut, latitude=latitude, nEVP=257);
[F, G, h, omega] = im.ModesAtWavenumber(2*pi/1000);
```




## Topics
+ Create and initialize modes
  + [`InternalModesSpectral`](/internal-modes/classes/numerical-solvers/internalmodesspectral/internalmodesspectral.html) Initialize the spectral solver on depth coordinates.
  + [`nEVP`](/internal-modes/classes/numerical-solvers/internalmodesspectral/nevp.html) Number of collocation points used for the generalized EVP.
+ Compute modes
  + [`BottomModesAtWavenumber`](/internal-modes/classes/numerical-solvers/internalmodesspectral/bottommodesatwavenumber.html) Return the bottom SQG mode at fixed horizontal wavenumber.
  + [`BoundaryModesAtWavenumber`](/internal-modes/classes/numerical-solvers/internalmodesspectral/boundarymodesatwavenumber.html) Return either the surface or bottom boundary mode at fixed wavenumber.
  + [`EigenmatricesForFrequency`](/internal-modes/classes/numerical-solvers/internalmodesspectral/eigenmatricesforfrequency.html) Assemble the fixed-`\omega` generalized EVP on the spectral grid.
  + [`EigenmatricesForWavenumber`](/internal-modes/classes/numerical-solvers/internalmodesspectral/eigenmatricesforwavenumber.html) Assemble the fixed-`K` generalized EVP on the spectral grid.
  + [`SurfaceModesAtWavenumber`](/internal-modes/classes/numerical-solvers/internalmodesspectral/surfacemodesatwavenumber.html) Return the surface SQG mode at fixed horizontal wavenumber.


## Developer Topics
These items document internal implementation details and are not part of the primary public API.
+ Configure normalization and boundaries
  + [`ApplyBoundaryConditions`](/internal-modes/classes/numerical-solvers/internalmodesspectral/applyboundaryconditions.html) Apply the active surface and bottom conditions to an EVP pair.
+ Developer topics
  + [`ChebyshevDifferentiationMatrix`](/internal-modes/classes/numerical-solvers/internalmodesspectral/chebyshevdifferentiationmatrix.html) Return the matrix that differentiates Chebyshev coefficients.
  + [`ChebyshevInterpolationDerivative`](/internal-modes/classes/numerical-solvers/internalmodesspectral/chebyshevinterpolationderivative.html) Return the Lobatto-grid differentiation matrix for interpolation values.
  + [`ChebyshevPolynomialsOnGrid`](/internal-modes/classes/numerical-solvers/internalmodesspectral/chebyshevpolynomialsongrid.html) Evaluate Chebyshev polynomials and derivatives on an arbitrary grid.
  + [`ChebyshevTransformForGrid`](/internal-modes/classes/numerical-solvers/internalmodesspectral/chebyshevtransformforgrid.html) Build a spectral interpolation map from a Lobatto grid to an output grid.
  + [`CheckIfReasonablyMonotonic`](/internal-modes/classes/numerical-solvers/internalmodesspectral/checkifreasonablymonotonic.html) Test whether a gridded density profile is monotonic enough for stretched coordinates.
  + [`Diff1_xCheb`](/internal-modes/classes/numerical-solvers/internalmodesspectral/diff1_xcheb.html) First-derivative operator in Chebyshev coefficient space.
  + [`DifferentiateChebyshevVector`](/internal-modes/classes/numerical-solvers/internalmodesspectral/differentiatechebyshevvector.html) Differentiate a Chebyshev series in coefficient space.
  + [`EigenmatricesForGeostrophicFModes`](/internal-modes/classes/numerical-solvers/internalmodesspectral/eigenmatricesforgeostrophicfmodes.html) Assemble the geostrophic interior `F`-mode EVP.
  + [`EigenmatricesForGeostrophicGModes`](/internal-modes/classes/numerical-solvers/internalmodesspectral/eigenmatricesforgeostrophicgmodes.html) Assemble the geostrophic interior `G`-mode EVP.
  + [`EigenmatricesForGeostrophicRigidLidGModes`](/internal-modes/classes/numerical-solvers/internalmodesspectral/eigenmatricesforgeostrophicrigidlidgmodes.html) Assemble the rigid-lid geostrophic EVP with displaced boundaries.
  + [`EigenmatricesForMDAModes`](/internal-modes/classes/numerical-solvers/internalmodesspectral/eigenmatricesformdamodes.html) Assemble the MDA diagnostic EVP.
  + [`EigenmatricesForRigidLidGModes`](/internal-modes/classes/numerical-solvers/internalmodesspectral/eigenmatricesforrigidlidgmodes.html) Assemble the rigid-lid geostrophic `G`-mode EVP.
  + [`EigenmatricesForSmithVannesteModes`](/internal-modes/classes/numerical-solvers/internalmodesspectral/eigenmatricesforsmithvannestemodes.html) Assemble the Smith-Vanneste SQG-like EVP used by a helper workflow.
  + [`FFromVCheb`](/internal-modes/classes/numerical-solvers/internalmodesspectral/ffromvcheb.html) Map Chebyshev coefficients to `F_j` on the active internal grid.
  + [`FNorm`](/internal-modes/classes/numerical-solvers/internalmodesspectral/fnorm.html) `\omega`-constant normalization functional.
  + [`FOutFromVCheb`](/internal-modes/classes/numerical-solvers/internalmodesspectral/foutfromvcheb.html) Map Chebyshev coefficients to `F_j(z)` on the public output grid.
  + [`FindRootsFromChebyshevVector`](/internal-modes/classes/numerical-solvers/internalmodesspectral/findrootsfromchebyshevvector.html) Find physical-domain roots of a Chebyshev series.
  + [`FindSmallestChebyshevGridWithNoGaps`](/internal-modes/classes/numerical-solvers/internalmodesspectral/findsmallestchebyshevgridwithnogaps.html) Find the coarsest Lobatto grid that covers a target grid without gaps.
  + [`FindTurningPointBoundariesAtFrequency`](/internal-modes/classes/numerical-solvers/internalmodesspectral/findturningpointboundariesatfrequency.html) Find turning points and region signs for `N2(z)-\omega^2`.
  + [`GFromVCheb`](/internal-modes/classes/numerical-solvers/internalmodesspectral/gfromvcheb.html) Map Chebyshev coefficients to `G_j` on the active internal grid.
  + [`GNorm`](/internal-modes/classes/numerical-solvers/internalmodesspectral/gnorm.html) `K`-constant normalization functional.
  + [`GOutFromVCheb`](/internal-modes/classes/numerical-solvers/internalmodesspectral/goutfromvcheb.html) Map Chebyshev coefficients to `G_j(z)` on the public output grid.
  + [`GaussQuadraturePointsForEigenmatrices`](/internal-modes/classes/numerical-solvers/internalmodesspectral/gaussquadraturepointsforeigenmatrices.html) Return quadrature points inferred from a generalized EVP.
  + [`GaussQuadraturePointsForMDAModes`](/internal-modes/classes/numerical-solvers/internalmodesspectral/gaussquadraturepointsformdamodes.html) Return quadrature points tailored to the diagnostic MDA modes.
  + [`GaussQuadraturePointsForModesAtFrequency`](/internal-modes/classes/numerical-solvers/internalmodesspectral/gaussquadraturepointsformodesatfrequency.html) Return quadrature points tailored to fixed-`\omega` modes.
  + [`GaussQuadraturePointsForModesAtWavenumber`](/internal-modes/classes/numerical-solvers/internalmodesspectral/gaussquadraturepointsformodesatwavenumber.html) Return quadrature points tailored to fixed-`K` modes.
  + [`GeostrophicFModesAtWavenumber`](/internal-modes/classes/numerical-solvers/internalmodesspectral/geostrophicfmodesatwavenumber.html) Return the geostrophic interior `F`-modes at fixed horizontal wavenumber.
  + [`GeostrophicModesAtWavenumber`](/internal-modes/classes/numerical-solvers/internalmodesspectral/geostrophicmodesatwavenumber.html) Return the geostrophic interior `G`-modes at fixed horizontal wavenumber.
  + [`GeostrophicNorm`](/internal-modes/classes/numerical-solvers/internalmodesspectral/geostrophicnorm.html) Geostrophic normalization functional.
  + [`GeostrophicRigidLidModes`](/internal-modes/classes/numerical-solvers/internalmodesspectral/geostrophicrigidlidmodes.html) Return rigid-lid geostrophic modes with displaced boundaries.
  + [`GeostrophicSmithVannesteModesAtWavenumber`](/internal-modes/classes/numerical-solvers/internalmodesspectral/geostrophicsmithvannestemodesatwavenumber.html) Return Smith-Vanneste modes at fixed horizontal wavenumber.
  + [`Int_xCheb`](/internal-modes/classes/numerical-solvers/internalmodesspectral/int_xcheb.html) Exact integral weights for Chebyshev coefficients on the active domain.
  + [`IntegrateChebyshevVector`](/internal-modes/classes/numerical-solvers/internalmodesspectral/integratechebyshevvector.html) Integrate a Chebyshev series in coefficient space.
  + [`IntegrateChebyshevVectorWithLimits`](/internal-modes/classes/numerical-solvers/internalmodesspectral/integratechebyshevvectorwithlimits.html) Integrate a Chebyshev series between two physical limits.
  + [`IsChebyshevGrid`](/internal-modes/classes/numerical-solvers/internalmodesspectral/ischebyshevgrid.html) Test whether a grid is a Chebyshev-Lobatto grid up to tolerance.
  + [`Lx`](/internal-modes/classes/numerical-solvers/internalmodesspectral/lx.html) Total span of the active spectral coordinate.
  + [`MDAModes`](/internal-modes/classes/numerical-solvers/internalmodesspectral/mdamodes.html) Return the diagnostic MDA modes.
  + [`N2_function`](/internal-modes/classes/numerical-solvers/internalmodesspectral/n2_function.html) Buoyancy-frequency profile represented as a function of depth.
  + [`N2_xLobatto`](/internal-modes/classes/numerical-solvers/internalmodesspectral/n2_xlobatto.html) Buoyancy frequency squared sampled on `xLobatto`.
  + [`ProjectOntoChebyshevPolynomialsWithTolerance`](/internal-modes/classes/numerical-solvers/internalmodesspectral/projectontochebyshevpolynomialswithtolerance.html) Project a profile onto Chebyshev coefficients to a requested tolerance.
  + [`T_xCheb_zOut`](/internal-modes/classes/numerical-solvers/internalmodesspectral/t_xcheb_zout.html) Transform from Chebyshev coefficients to the public output grid.
  + [`T_xLobatto`](/internal-modes/classes/numerical-solvers/internalmodesspectral/t_xlobatto.html) Chebyshev basis evaluated on the Lobatto grid.
  + [`Tx_xLobatto`](/internal-modes/classes/numerical-solvers/internalmodesspectral/tx_xlobatto.html) First derivatives of the Chebyshev basis on the Lobatto grid.
  + [`Txx_xLobatto`](/internal-modes/classes/numerical-solvers/internalmodesspectral/txx_xlobatto.html) Second derivatives of the Chebyshev basis on the Lobatto grid.
  + [`ValueOfFunctionAtPointOnGrid`](/internal-modes/classes/numerical-solvers/internalmodesspectral/valueoffunctionatpointongrid.html) Evaluate a Chebyshev series at arbitrary points in its physical domain.
  + [`fInverseBisection`](/internal-modes/classes/numerical-solvers/internalmodesspectral/finversebisection.html) Invert a monotonic function by vectorized bisection.
  + [`fct`](/internal-modes/classes/numerical-solvers/internalmodesspectral/fct.html) Compute the fast Chebyshev transform of values on a Lobatto grid.
  + [`hFromLambda`](/internal-modes/classes/numerical-solvers/internalmodesspectral/hfromlambda.html) Set on initialization by the subclass, these transformations are
  + [`ifct`](/internal-modes/classes/numerical-solvers/internalmodesspectral/ifct.html) Compute the inverse fast Chebyshev transform.
  + [`rho_function`](/internal-modes/classes/numerical-solvers/internalmodesspectral/rho_function.html) Density profile represented as a spline or functional object.
  + [`standardChop`](/internal-modes/classes/numerical-solvers/internalmodesspectral/standardchop.html) Choose a truncation index for a Chebyshev coefficient sequence.
  + [`xDomain`](/internal-modes/classes/numerical-solvers/internalmodesspectral/xdomain.html) Bounds of the active spectral coordinate.
  + [`xLobatto`](/internal-modes/classes/numerical-solvers/internalmodesspectral/xlobatto.html) The 'x' refers to the stretched coordinate being used.
  + [`xMax`](/internal-modes/classes/numerical-solvers/internalmodesspectral/xmax.html) Maximum value of the active spectral coordinate.
  + [`xMin`](/internal-modes/classes/numerical-solvers/internalmodesspectral/xmin.html) Minimum value of the active spectral coordinate.
  + [`xOut`](/internal-modes/classes/numerical-solvers/internalmodesspectral/xout.html) Output locations mapped into the active spectral coordinate.
  + [`x_function`](/internal-modes/classes/numerical-solvers/internalmodesspectral/x_function.html) Active stretched-coordinate map `x(z)` used by the solver.
  + [`z_xLobatto`](/internal-modes/classes/numerical-solvers/internalmodesspectral/z_xlobatto.html) Physical depths associated with `xLobatto`.


---