# Continuous basis contracts

`IMModes` resolves a physical mode request into an immutable `IMModeBasis`. The basis exposes three related value records: `scientificSpecification`, `constructionRecord` and `capabilities`. A collection copies these records from each common basis. Creating a record or selecting existing columns performs no additional EVP solve, grid design or quadrature fit.

## Scientific definition and construction

`scientificSpecification` records the family and background, the requested scientific parameters and their resolved meaning. Its `resolved` fields include:

| Field | Meaning |
| --- | --- |
| `evp` | Canonical source formulation, coefficient functions, parameters and boundary coefficients, when the family has an EVP |
| `pairings` | Defining physical or label functionals; these describe the pairing, rather than a sampled or normalized Gram matrix |
| `energyBoundary`, `energyBoundaryRole` | Generalized-energy coefficients and their role, separately from `background.surfaceCondition` |
| `intervals`, `intervalKind`, `intervalBoundaryTypes`, `endpointDerivativeTrace` | Declared neutral/critical support and its endpoint interpretation |
| `profileRelation`, `profileBoundarySemantics` | Common F/G/APV convention and source boundary semantics |
| `normalization` | Name, positive divisors `scales` and separate `normSignatures` |
| `horizontal` | Fixed, per-column or independent horizontal identity |
| `columnLabels`, `modeInfo`, `sourceColumns` | Current retained order and one-to-one source indices when applicable; mixed zero-APV columns instead use the explicit source-coordinate map |
| `includesAnalyticalBarotropic` | Whether the current hydrostatic-F inventory contains the analytical barotropic column |

Normalized profiles equal the unscaled common-family profiles, after applying `profileRelation`, divided by `normalization.scales`, column by column. For waves this includes replacing the source F convention with F=DG-prime before scaling. Every field in a column has the same divisor. For zero-APV bases, `sourceCoordinateMap` also composes rotations and selections from `canonicalEndpointLabels` into current normalized columns. Signed metric signatures remain separate. `basis.metadata` is reserved for caller extensions and is not the authority for these records.

The physical surface convention and generalized-energy coefficients are independent inputs. APV and mean-density modes use the coefficients to define their pairing and endpoint coordinates. The canonical zero-APV response profiles use their active endpoint identities. For neutral DD labels the optional coefficients are provenance only: G vanishes at both ends of every interval, including physical boundaries. Interval derivative values use the inside trace; zero extension does not imply a globally classical derivative across a support endpoint.

`constructionRecord.requested.solver` preserves the backend request. `resolved` identifies the actual backend/source class and available solver settings, including the configured native representation when a source retains it. Analytical barotropic and endpoint-empty bypasses explicitly report skipped construction. Empty numerical selections can discard their native representation; `nativeRepresentationAvailable=false` reports that absence without inventing nodes or reconfiguring the solver.

`constructionRecord.selection` keeps the original ladder request and backend selection diagnostics, along with the current `activeSourceColumns`, `columnLabels` and composed `initialColumnIndices` when the mapping is one-to-one. `sourceEigenpairIndices` identifies the backend eigenpairs; `activeSourceColumns` addresses the extended inventory, which also includes the analytical barotropic column when present. Dense zero-APV combinations use an explicit coordinate map rather than fictitious single source indices. `selectColumns` updates the active identities, spectra and normalization divisors while preserving the original construction request and diagnostics. Requested-empty and subsequently selected-empty bases retain their declared operations.

`constructionRecord.evidence.solveAccuracy` and `referenceConvergence` are explicitly `"unverified"` unless a separate qualification supplies evidence. A correct returned count, a native resolution and an analytical backend name do not substitute for that evidence. Source diagnostics are retained separately.

## Horizontal identity

| `horizontalKind` | Families | `k` | `kSquared` |
| --- | --- | --- | --- |
| `"fixedWavenumber"` | Fixed-k waves, inertial auxiliary modes, generalized-enstrophy modes, zero-APV responses and zero-depth waves | One nonnegative solved k | Empty |
| `"perColumn"` | Fixed-frequency waves | Empty | One signed horizontal root per column |
| `"independent"` | Hydrostatic F/G, APV, mean-density anomaly and neutral labels | Empty | Empty |
| `"unspecified"` | Developer wrappers without an explicit horizontal declaration | Empty | Empty |

The genuine inertial construction has `k=0`. An independent basis does not acquire a fictitious zero wavenumber.

For fixed-frequency modes, the source eigenvalue is the inverse equivalent depth, $$\lambda_j=1/h_j$$. Horizontal dispersion is

$$k_j^2=\frac{(\omega-|f_0|)(\omega+|f_0|)}{g}\lambda_j.$$

Positive roots propagate; negative roots are evanescent. The sign of h alone is insufficient, because the frequency factor changes sign below the inertial frequency. Keeping signed `kSquared` also preserves a zero root without taking a reciprocal of h. The inertial-frequency degeneracy uses the separate inertial factory.

An `IMBasisCollection` accepts caller-supplied `kappa` for explicitly independent bases, so the same vertical structure can serve positive horizontal wavenumbers. Fixed-k pages must match their solved k exactly. Per-column roots and unspecified wrappers reject assignment of a scalar page kappa. With omitted `kappa`, the collection infers a row only when every mapped page has a known scalar wavenumber; otherwise `collection.kappa` remains empty. Each basis's horizontal identity is still available in `collection.metadata`.

```matlab
background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
apv = IMModes.geostrophicAPVModes(background, ...
    energyBoundary=IMGeneralizedEnergyBoundary(g0=Inf,gd=Inf), ...
    solver=IMSolverAnalytical(),nPositiveModes=4);
pages = IMBasisCollection({apv},kappa=[1e-3 2e-3 1e-3],basisIndex=[1 1 1]);
```

## Shared capability authority

`basis.capabilities` is a structure, replacing the former variable-name row in this prerelease API. Its canonical `variables` and parallel `derivativeOrders` cells describe evaluable channels. Aliases `Fz`, `Gz` and `Gzz` map onto a canonical variable and derivative offset; derivative limits still apply after resolving an alias. Supported derivatives are physical vertical derivatives, with interval trace conventions explicit. Scalar evaluation, collection evaluation, chunked evaluation and product factors consume this same record for common bases. Unsupported requests fail before profile evaluation, including on an empty inventory.

```matlab
description = apv.capabilities;
description.variables
description.derivativeOrders
Fz = apv.evaluate([-900;-500;-100],"Fz");
```

`spectrumNames` distinguishes available spectra from absent ones, including empty spectra. `normalizationNames` lists structurally implemented rules. A retained null or nonfinite column can still make a particular normalization inadmissible; this is reported rather than dropping that column. A scientific pairing can exist before a particular backend implements its normalization operation. Direct intrinsic wrappers declare their implemented variables/derivatives explicitly and do not gain F-prime merely by naming an F profile.

`basis.withNormalization(name)` returns another immutable basis using the same continuous representation and selected columns, with updated resolved divisors and signatures. It preserves the original construction requests and evidence. See [changing normalization](BasisNormalization.md) for admissibility, selection invariance and the corresponding changes to coefficients, physical duals and recovery operators.

Zero-APV bases also expose [physical coefficient forms and rotations](BasisRotations.md). Form getters use current normalized/selected coordinates. Canonical endpoint labels identify source rows; rotated column labels identify their linear combinations. Collection metadata preserves both sets and their map. Its `columnKind="endpoint"` identifies the finite endpoint sector and its fixed-count assessment semantics, rather than asserting that each rotated column belongs to one endpoint.

`projectionOperations` records the current conditional physical adapters: rigid-lid hydrostatic-F provides sampled Galerkin F/G projection, and neutral labels provide the prescribed continuous APV dual. Requirements remain explicit. Other families can expose a defining pairing without yet having a common physical projection adapter.

`sampledRecovery` describes synthesis channels and the requirements for preparing numerical recovery. It does not claim rank, noise control, completeness or recovery success for any observation grid. Fitting evidence, observation/noise factors, activity declarations and measured rank belong to the prepared transform.

## Grid and recovery ownership

The scientific EVP, solver coordinate, observation grid and numerical recovery are separate choices. A provider-designed root grid is one observation geometry; an externally imposed model grid is another. Both use the same continuous identities and can feed the [shared recovery engine](SharedRecovery.md). A diagnostic application may evaluate the required channels on its saved-output grid and pass its assembled observation matrices and factors to `IMDiscreteTransform.fromObservationMatrices`. Boundary traces, derivatives, missing channels and particle labels must be represented by the actual observation model.

This contract does not choose an observation grid, qualify an inverse or define WVM polarization. Additional common physical projection adapters, rotations outside the supported zero-APV sector, model-specific channel assembly and WVM adoption remain separate increments. See [mode creation](ModeCreationAPI.md), [vertical grids](VerticalGrids.md) and [hydrostatic root transforms](HydrostaticRootTransform.md).
