# Zero-APV forms and basis rotations

`IMModes.geostrophicZeroAPVModesAtWavenumber` returns a common continuous basis with one canonical response for each active endpoint. The basis also exposes physical forms and immutable rotations of its retained subspace. These operations reuse the continuous representation; they require no new boundary-value solve, observation grid or weight fit.

```matlab
background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
boundary = IMGeneralizedEnergyBoundary(g0=-.035,gd=.01);
responses = IMModes.geostrophicZeroAPVModesAtWavenumber(background,1e-3,energyBoundary=boundary);
Hg = responses.generalizedEnergyMatrix();
depthModes = responses.rotateBoundaryDepth();
rmsModes = depthModes.withNormalization("depthRMSF");
```

## Physical forms in current coordinates

| API | Restricted physical form |
| --- | --- |
| `energyMatrix` | Physical energy H |
| `surfaceBuoyancyMatrix` | Scaled surface-buoyancy variance B0 |
| `bottomBuoyancyMatrix` | Scaled bottom-buoyancy variance Bd |
| `endpointResponseMetric` | Sum of squared active endpoint responses RB |
| `generalizedEnergyMatrix(g0=...,gd=...)` | Generalized energy Hg = H + g0 B0 + gd Bd |

Every matrix has `nColumns` rows and columns and uses the basis's current normalized, selected or rotated coordinates. The forms retain the existing lower-level zero-APV state convention, including its factors of f0, g and k. These are physical coefficient forms, rather than raw depth integrals of a single profile. Physical energy of the boundary-buoyancy component is a self-pairing; it need not be an additive reservoir when combined with PV modes, because physical-energy cross terms can remain.

For current profile rows Fs, Fb, Gs and Gb at the surface and bottom, let bs = Gs - Fs for a free surface and bs = Gs for a rigid lid, with bd = Gb. The restricted forms are

$$H=\frac{f_0^2}{2g k^4}\operatorname{sym}(-F_s^*b_s+F_b^*b_d),\qquad B_0=\frac{f_0^2}{2g^2 k^4}b_s^*b_s,\qquad B_d=\frac{f_0^2}{2g^2 k^4}b_d^*b_d,$$

where `sym` takes the symmetric part. The response metric is the sum of the active endpoint response outer products, and an absent endpoint's buoyancy form is exactly zero. These are the manuscript's boundary-buoyancy restricted state forms, including its distinction between total displacement and the free-surface buoyancy response.

Omitted `g0` and `gd` use the retained scientific energy-boundary coefficients. Explicit overrides change the requested form or rotation and are recorded with the rotation. Terms for endpoints already absent from `canonicalEndpointLabels` are omitted, including an infinite inactive-endpoint coefficient. An override cannot deactivate an endpoint retained by the source; an infinite coefficient there is rejected. Signed finite coefficients remain signed, and a null generalized-energy direction cannot be assigned unit generalized energy.

Zero-APV generalized-energy normalization rejects an unresolved self-pairing at relative tolerance `1e-12`. It compares the signed magnitude with both the physical majorant H + |g0| B0 + |gd| Bd and a per-column cancellation bound from the entrywise absolute canonical forms and source directions. The latter detects roundoff in a physically null neutral direction even when a rounded physical-majorant self-pairing is tiny and positive. There is no dimensional absolute floor. These tests reject an unresolved norm; they do not turn the signed form into a positive pairing or remove the direction from the basis. A suitable profile normalization remains available.

## Canonical endpoints and current columns

`canonicalEndpointLabels` gives the source-row identities, such as `["surface","bottom"]`. `columnLabels` identifies current basis columns. They coincide for the canonical responses but have different meanings after mixing columns.

The rectangular `sourceCoordinateMap` M is authoritative:

$$A_{\mathrm{current}} = A_{\mathrm{canonical}} M.$$

This equation holds for F, G, their supported physical derivatives and endpoint traces. Thus canonical endpoint coefficients equal M times current coefficients, and every form changes by congruence:

$$H_{\mathrm{current}} = M^* H_{\mathrm{canonical}} M.$$

Selecting columns restricts M's columns. Normalization rescales them by positive divisors. A rotation composes M on the right with a current-coordinate transformation. A mixed column has no single source index or endpoint identity. Current coordinate records preserve the map rather than inventing a one-to-one source-column assignment; original construction requests and backend evidence remain available.

`scientificSpecification.resolved.sourceCoordinateMap` records current normalized coordinates. `constructionRecord.selection.sourceDirections` records the fixed directions before subsequent profile divisors and therefore stays unchanged by renormalization. New rotated identities clear `initialColumnIndices`; their source directions provide the appropriate provenance.

The developer constructor can wrap an already rotated lower-level source and retains its canonical map. It infers an omitted source wavenumber and rejects contradictory fixed k, g, f0, domain or surface convention. The supplied background must also describe the same stratification as the source; equivalence of arbitrary coefficient functions is a developer precondition rather than a numerical certificate. The `IMModes` factories construct the source from the declared background directly.

## Supported rotations

All pencils below act in the current retained coordinates. A selected subspace stays selected; rotation never restores discarded source directions.

| Method | Left form L | Right form R | Rotation normalization |
| --- | --- | --- | --- |
| `rotateBoundaryDepth()` | Hg | RB | Unit endpoint-response metric |
| `rotateSurfaceBuoyancy()` | g B0 | Hg | Unit signed generalized energy |
| `rotateWithPencil(name=...,leftMatrix=L,rightMatrix=R)` | Caller-supplied current-coordinate L | Caller-supplied current-coordinate R | Unit signed right metric |

The real symmetric pencil satisfies

$$L T = R T \operatorname{diag}(\lambda),\qquad T^* R T = \operatorname{diag}(\sigma),\qquad \sigma_j\in\{-1,1\}.$$

The returned basis has synthesis `Anew=Aold*T`. Boundary-depth eigenvalues are ordered descending and define signed boundary depths $$h_0=2k^2\lambda$$. Surface-buoyancy directions with nonzero eigenvalue come first and require an active surface endpoint. Custom matrices must have `nColumns` rows and columns; their meaning follows the common form getters. The older lower-level custom-pencil API instead takes canonical boundary coordinates, so copying its matrices after changing common coordinates requires the corresponding congruence.

`rotationMetadata` records the name, eigenvalues, metric signatures, normalization convention, residuals, boundary depths when defined and the transformation from the parent coordinates. These describe the rotation that was performed. They remain distinct from the current profile normalization after a subsequent `withNormalization` call. Zero-APV modes have no differential-EVP eigenvalues or equivalent depths: `basis.eigenvalues` and `basis.h` remain unavailable.

The immediate rotation normalization is `"rotationRightMetric"`; its exact physical meaning is in `rotationMetadata.normalizationConvention`. A unit right-metric vector is generally neither a canonical unit endpoint response nor a depth-RMS profile. Later positive profile normalization retains the coordinate map and the original rotation evidence. Call `withNormalization("rotationRightMetric")` to restore that rotation's retained right-metric normalization. Signed metrics are never replaced by their absolute-value forms to manufacture an admissible pencil.

The supported pencils require a complete, finite, real, numerically independent eigenbasis with nonnull right-metric norms. Singular, defective, nonreal or unresolved metric cases fail explicitly. A zero left eigenvalue is allowed when the right metric remains admissible; its generalized-energy norm may still be zero. Empty inventories preserve empty forms and coordinate maps.

## Prepared recovery and qualification

For an invertible change T on the same active span, physical duals B, prepared recovery R and coefficient scaling C obey

$$B_{\mathrm{new}}=T^{-1}B,\qquad R_{\mathrm{new}}=T^{-1}R,\qquad C_{\mathrm{new}}=T^{-1}C.$$

The observation and fitting factors, recovery policy, active span and retained/extended interpretation must retain their physical meanings. These relations also apply to inputs outside the represented span. A restriction to fewer columns changes the represented subspace and is not an invertible coordinate equivalence. General rotations across active and inactive sectors are outside this contract.

This API supplies continuous forms and coordinate changes to the [shared recovery engine](SharedRecovery.md). It does not add a separate recovery engine or claim a new family/source projection adapter. See [verification](BasisRotationVerification.md) for tested cases and limitations; broader adapters remain a separate increment.
