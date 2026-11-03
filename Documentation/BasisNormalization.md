# Changing common-basis normalization

`basis.withNormalization(name)` returns an independent basis value using the existing continuous representation and retained columns. It supports the rules listed in `basis.capabilities.normalizationNames`, subject to each retained column having an admissible magnitude.

```matlab
background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
problem = IMInternalModes.waveModesAtWavenumber(background,k=1e-3);
waves = IMSolverAnalytical().solveEVP(problem,nModes=9);
rmsWaves = waves.withNormalization("rmsF");
```

`waves` keeps its original normalization. `rmsWaves` preserves the same physical background, source representation, scientific labels, spectra, horizontal identity and selected column order. Every related field and physical derivative in a column has the same new positive divisor. Profile relations such as wave $$F=D G'$$ are preserved.

Normalization factors are computed lazily when profiles, forms or factors are requested; merely changing a supported default need not evaluate them. It reuses the existing modes and performs no physical EVP solve, observation-grid construction or quadrature-weight fit. Requesting the current normalization returns the existing value directly. Previously prepared transforms keep their captured basis. When preparing a new adapter on the same grid, preserve the original recovery policy, retained/extended spans and activity to obtain an equivalent recovery map.

## Meaning and admissibility

Normalization names retain their [mode-creation definitions](ModeCreationAPI.md#profile-convention-and-lazy-normalization). Names describe implemented rules, and a particular retained inventory must also be admissible. For example, the analytical barotropic column has G=0, so an inventory containing it cannot use a G-based unit normalization. Select an admissible subspace explicitly before requesting that rule. No column is silently removed.

All built-in normalizations keep the profiles dimensionless. `rmsF`, `rmsG` and generic scalar `rms` use a full-depth average, without boundary contributions. Their names refer to profiles: F is not always horizontal velocity, and G is not always a displacement profile.

Every canonical name is reserved across families. Register a distinct custom name; that custom rule can be replaced without affecting the source value. V2 accepts text names only and rejects legacy aliases, `Normalization` enums, `innerProduct*` variations and `generalizedEnergy`. Energy matrices and `orthogonalizeGeneralizedEnergy` retain their physical units and signed meaning independently of profile normalization.

`meanDensityAnomalyEnergy` preserves the manuscript's dimensionless signed G form with magnitude one, so its physical generalized-energy self-pairing is g times its signature divided by two. The new name preserves the earlier MDA scale. Null or nonfinite magnitudes fail explicitly.

Unsupported rules are rejected before profile evaluation, including on empty inventories. Supported empty inventories retain empty scales and signatures, with their capability and scientific identity records intact.

Scientific records identify the active normalization without forcing its integrals. For endpoint bases the current `boundaryResponseMatrix` includes the active column scales; orthogonalization evidence and intrinsic `boundaryDepth` retain their original meaning. Original scientific requests and backend construction evidence remain available. Selecting and reordering columns preserves these distinctions. Normalization integration retains the original source-inventory budget through selection, so selecting a small high-order band does not reduce the integration rule to that of a few low modes. This fixes a selected-harmonic aliasing error; it is not an independent convergence certificate for arbitrary supplied profiles.

For boundary-buoyancy modes, choose [canonical or generalized-energy orthogonal coordinates](BasisRotations.md) before reducing the inventory. Profile normalization is a separate operation. The default orthogonal basis uses the manuscript's unit endpoint-response geometry. A null generalized-energy direction remains a valid basis vector; an inverse-energy projection on that direction is unavailable. Explicit coefficients define that pairing independently of endpoint selection.

## Coefficient coordinates and prepared operators

Let the original and new positive divisors be `sOld` and `sNew`, in the same scientific column order, and define the representable coordinate change with componentwise division

$$T=\operatorname{diag}(s_\mathrm{old}/s_\mathrm{new}).$$

Every channel's synthesis matrix changes as $$A_\mathrm{new}=A_\mathrm{old}T$$. A fixed physical field consequently has coefficients $$a_\mathrm{new}=T^{-1}a_\mathrm{old}$$. Coordinate-dependent physical duals B, numerical recovery maps R and physical coefficient forms H satisfy

$$B_\mathrm{new}=T^{-1}B_\mathrm{old},\qquad R_\mathrm{new}=T^{-1}R_\mathrm{old},\qquad H_\mathrm{new}=T^*H_\mathrm{old}T.$$

For caller-prepared [shared recovery](SharedRecovery.md), also transform the coefficient-coordinate map as $$C_\mathrm{new}=T^{-1}C_\mathrm{old}$$. The physical observation and fitting factors retain their meaning and values. Recovery covariance compares the same policy, retained/extended spans and activity; grid points and weights alone do not specify these choices. Consistent scaling preserves the physical matrix problem, with numerical preparation still subject to the engine's representability and rank checks. Comparing coefficient numbers or raw dimensional condition numbers alone does not establish physical equivalence.

In particular, the root factory defaults to `"truncatedCollocation"`, using recovery against the full root span. Calling `normalizedBasis.discreteTransform(grid)` alone prepares retained-span Galerkin recovery. To preserve the original root recovery, carry its extended-span specification through the shared engine and pass the corresponding prepared F/G recoveries to the adapter. This basis operation does not mutate or automatically convert an existing transform.

These relations apply within the supported active coordinates. Positive diagonal normalization preserves activity, including the inactive hydrostatic G-barotropic column. Scientific inactivity remains explicit; a sampled-zero active mode still fails recovery. General rotations that mix active and inactive sectors need a separate contract.

Physical dual application remains distinct from fitted recovery. Verification includes source inputs outside the retained modal span, for which these operations can differ even though each changes coordinates consistently. For the prescribed neutral APV dual, its identity duality target remains identity; its physical APV majorant Gram changes by congruence.

## Scope and verification

Normalization adds no recovery policy. Supported scalar physical adapters are described in [scalar projection](CommonScalarProjection.md); observation-grid qualification and complete WVM state/source assembly remain separate responsibilities.

The [historical normalization verification](BasisNormalizationVerification.md) records the incoming implementation. The integration verification must additionally cover capable families, lazy rules, retained-only admissibility and cache independence.
