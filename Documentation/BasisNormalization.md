# Changing common-basis normalization

`IMModeBasis.withNormalization(name)` returns a new immutable basis using the existing continuous representation and retained columns. It supports the rules listed in `basis.capabilities.normalizationNames`, subject to each retained column having an admissible magnitude.

```matlab
background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
waves = IMModes.waveModesAtWavenumber(background,1e-3,solver=IMSolverAnalytical(),nInternalModes=8);
rmsWaves = waves.withNormalization("depthRMSF");
```

`waves` keeps its original normalization. `rmsWaves` preserves the same physical background, source representation, scientific labels, spectra, horizontal identity and selected column order. Every related field and physical derivative in a column has the same new positive divisor. Profile relations such as wave $$F=D G'$$ are preserved.

Changing normalization can evaluate profiles and continuous normalization integrals. It reuses the existing modes and performs no physical EVP solve, observation-grid construction or quadrature-weight fit. Requesting the current normalization returns the existing value directly. Previously prepared transforms keep their captured basis. When preparing a new adapter on the same grid, preserve the original recovery policy, retained/extended spans and activity to obtain an equivalent recovery map.

## Meaning and admissibility

Normalization names retain their [mode-creation definitions](ModeCreationAPI.md#profile-convention-and-normalization). Names describe implemented rules, and a particular retained inventory must also be admissible. For example, the analytical barotropic column has G=0, so an inventory containing it cannot use a G-based unit normalization. Select an admissible subspace explicitly before requesting that rule. No column is silently removed.

Generalized-energy normalization uses the square root of the absolute signed self-pairing as its positive divisor and stores the sign separately in `normSignatures`. The underlying physical pairing retains its sign. Switching to a positive profile normalization changes the normalization signatures as appropriate; it does not replace the scientific signed pairing with a positive form. Null or nonfinite magnitudes are rejected explicitly.

Unsupported rules are rejected before profile evaluation, including on empty inventories. Supported empty inventories retain empty scales and signatures, with their capability and scientific identity records intact.

`scientificSpecification.resolved.normalization` records the new name, divisors and signatures. The original `scientificSpecification.requested` and `constructionRecord` remain the original construction evidence. Selecting and reordering columns preserves these distinctions. Normalization integration retains the original source-inventory budget through selection, so selecting a small high-order band does not reduce the integration rule to that of a few low modes. This fixes a selected-harmonic aliasing error; it is not an independent convergence certificate for arbitrary supplied profiles.

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

This increment provides the immutable basis operation and verifies the existing physical and recovery adapters under it. It adds no new recovery policy or family/source adapter. Zero-APV forms and rotations remain #49; common family/source adapter completion remains #48. Observation-grid qualification and WVM adoption remain separate work.

See [normalization verification](BasisNormalizationVerification.md) for numerical controls, independent review and the source-inventory integration correction.
