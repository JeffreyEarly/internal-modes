---
layout: default
title: CheckIfReasonablyMonotonic
parent: InternalModesSpectral
grand_parent: Classes
nav_order: 8
mathjax: true
---

#  CheckIfReasonablyMonotonic

Test whether a gridded density profile is monotonic enough for stretched coordinates.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [flag,dTotalVariation,rho_zCheb,rho_zLobatto,rhoz_zCheb,rhoz_zLobatto] = CheckIfReasonablyMonotonic(zLobatto,rho_zCheb,rho_zLobatto,rhoz_zCheb,rhoz_zLobatto)
```
## Parameters
+ `zLobatto`  Lobatto depth grid
+ `rho_zCheb`  Chebyshev coefficients of density
+ `rho_zLobatto`  density sampled on `zLobatto`
+ `rhoz_zCheb`  Chebyshev coefficients of $$\partial_z \rho$$
+ `rhoz_zLobatto`  $$\partial_z \rho$$ sampled on `zLobatto`

## Returns
+ `flag`  monotonicity status flag
+ `dTotalVariation`  fractional change in total variation after coercion
+ `rho_zCheb`  possibly corrected Chebyshev coefficients of density
+ `rho_zLobatto`  possibly corrected density on `zLobatto`
+ `rhoz_zCheb`  possibly corrected Chebyshev coefficients of $$\partial_z \rho$$
+ `rhoz_zLobatto`  possibly corrected $$\partial_z \rho$$ on `zLobatto`

## Discussion

                              We want to know if the density function is decreasing as z
  increases. If it's not, are the discrepencies small enough
  that we can just force them?

  flag is 0 if the function is not reasonably monotonic, 1 if
  it is, and 2 if we were able to to coerce it to be, without
  too much change
