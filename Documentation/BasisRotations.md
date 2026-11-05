# Boundary-buoyancy basis choices

Zero-APV geostrophic modes span the boundary-buoyancy states at each nonzero horizontal wavenumber. The manuscript uses two choices: canonical unit endpoint responses, and their generalized-energy orthogonal combinations. Both reuse the same continuous solutions.

```matlab
background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
problem = IMGeostrophicZeroAPVModes.atWavenumber(background,1e-3,endpoints=["surface","bottom"]);
responses = IMSolverSpectral(nEVP=64).solveGeostrophicZeroAPVModes(problem);
orthogonal = responses.orthogonalizeGeneralizedEnergy(g0=-.035,gd=.01);
hbb = orthogonal.boundaryDepth;
rmsModes = orthogonal.withNormalization("rmsF");
```

## Canonical endpoint coordinates

`endpoints` selects a possibly empty subset of surface and bottom. It is independent of the generalized-energy coefficients. Each canonical column has unit buoyancy response at its requested endpoint and zero response at the other endpoint. For a free surface the response is B(0)=G(0)-F(0); at the bottom B(-D)=G(-D). For a rigid lid the surface response is G(0).

`boundaryResponseMatrix` M maps current modal coefficients to canonical endpoint coefficients:

$$a_{\mathrm{boundary}}=M a_{\mathrm{current}},\qquad A_{\mathrm{current}}=A_{\mathrm{canonical}}M.$$

Its rows identify the originally requested canonical endpoints; columns identify the current modes; pages identify wavenumbers. Canonical M is the identity. Orthogonalization changes its directions, normalization rescales its columns, and selection restricts its columns. The same relation holds for F, G, B, derivatives and endpoint traces. No new boundary-value solve, observation grid or weight fit is needed.

Construct the complete endpoint inventory, orthogonalize if needed, then select and normalize. Orthogonalization requires the complete originally requested span. A proper selection, including an empty selection from a nonempty family, cannot be orthogonalized afterward. This prevents silently restoring a discarded direction. A single-endpoint family is complete. An originally empty family remains empty without an eigensolve.

## Generalized-energy orthogonal combinations

For explicitly supplied g0 and gd, form the real symmetric canonical generalized-energy matrix and diagonalize it in unit endpoint-response coordinates:

$$H_g C=C\operatorname{diag}\left(\frac{h_{\mathrm{bb}}}{2k^2}\right),\qquad C^T C=I.$$

`orthogonalizeGeneralizedEnergy(g0=...,gd=...)` returns the resulting combinations with M=C. The default scaling is the manuscript's endpoint-response scaling. `boundaryDepth` contains the signed manuscript quantities h_bb, rather than eigenvalues of the boundary-value differential equations. Canonical endpoint modes have no ordinary EVP equivalent depth or spectrum. Subsequent normalization changes physical self-pairings but preserves the intrinsic boundary-depth labels.

Negative, zero and repeated boundary depths are admissible. Generalized energy need not be positive definite, so its matrix is not used as a Cholesky metric. Numerical qualification uses relative symmetry, eigen-residual and coefficient-orthogonality checks. Repeated eigenvalues identify an invariant subspace, not unique individual vectors. Physical energy need not be diagonal in these coordinates.

Requested endpoints require finite signed coefficients. A positive infinite coefficient is permitted only at an absent endpoint, whose term is omitted before arithmetic. It never selects or removes endpoints. NaN and negative infinity are rejected.

## Physical forms in current coordinates

| API | Physical coefficient form |
| --- | --- |
| `energyMatrix` | Physical energy H |
| `surfaceBuoyancyMatrix` | Scaled surface-buoyancy variance B0 |
| `bottomBuoyancyMatrix` | Scaled bottom-buoyancy variance Bd |
| `generalizedEnergyMatrix(g0=...,gd=...)` | H_g = H + g0 B0 + gd Bd |

Every matrix uses the current selected and normalized coordinates. For current profile rows Fs, Fb, Gs and Gb, let bs=Gs-Fs for a free surface and bs=Gs for a rigid lid, with bd=Gb. The manuscript forms are

$$H=\frac{f_0^2}{2g k^4}\operatorname{sym}(-F_s^*b_s+F_b^*b_d),\qquad B_0=\frac{f_0^2}{2g^2 k^4}b_s^*b_s,\qquad B_d=\frac{f_0^2}{2g^2 k^4}b_d^*b_d.$$

An absent endpoint's buoyancy form is zero. Every form changes by congruence:

$$H_{\mathrm{current}}=M^*H_{\mathrm{canonical}}M.$$

 Physical energy can retain cross terms between PV and boundary-buoyancy states even when generalized energy makes those subspaces orthogonal.

Generalized-energy orthogonalization does not impose unit energy. The default remains `endpointResponse`, preserving dimensionless profiles and the manuscript's response geometry. A null generalized-energy direction remains in the basis; inverse-energy projection on it fails explicitly. Profile normalizations such as `rmsF` remain available when their own amplitudes are nonzero. The basis never removes a direction to make a norm succeed.

## Coefficients, source projection and recovery

For an invertible change T on the same active span, physical duals B and numerical recovery maps R transform as

$$B_{\mathrm{new}}=T^{-1}B,\qquad R_{\mathrm{new}}=T^{-1}R.$$

Keep observation factors, fitting policy and represented spans consistent when comparing prepared operators. Selecting fewer columns changes the physical subspace; prepare its projection in the retained coordinates rather than slicing a full-space physical dual.

For an ocean state, determine canonical endpoint coefficients from residual endpoint anomalies after subtracting the reconstructed PV contribution, as in the manuscript. Independent source observations remain independent inputs. `IMGeostrophicTransform` owns this composition and uses the response map to convert coordinates. A scalar profile alone does not define a complete state projector.

The common [normalization](BasisNormalization.md) and [shared recovery](SharedRecovery.md) contracts apply. The integration [plan](APISimplificationPlan.md) and its verification ledger record the API migration and numerical checks.
