# Continuous basis contracts

A scientific family problem and its solver return a capable numerical, analytical or specialized family basis. The basis exposes three related value records: `scientificSpecification`, `constructionRecord` and `capabilities`. A collection copies these records from each common basis. Creating a record or selecting existing columns performs no additional EVP solve, grid design or quadrature fit.

## Scientific definition and construction

`scientificSpecification` records the family and background, the requested scientific parameters and their resolved meaning. Its `resolved` fields include:

| Field | Meaning |
| --- | --- |
| `evp` | Canonical source formulation, coefficient functions, parameters and boundary coefficients, when the family has an EVP |
| `pairings` | Defining physical or label functionals; these describe the pairing, rather than a sampled or normalized Gram matrix |
| `energyBoundary`, `energyBoundaryRole` | Generalized-energy coefficients and their role, separately from `background.surfaceCondition` |
| `intervals`, `intervalKind`, `intervalBoundaryTypes`, `endpointDerivativeTrace` | Declared neutral/critical support and its endpoint interpretation |
| `profileRelation` | Common F/G/APV convention and source boundary semantics |
| `normalization` | Active name and lazy factor availability; factors and signatures are evaluated only when needed |
| `horizontal` | Fixed, per-column or independent horizontal identity |
| `columnLabels`, `modeInfo`, `sourceColumns` | Current retained order and one-to-one source indices when applicable; mixed zero-APV columns instead use the explicit source-coordinate map |
| `includesAnalyticalBarotropic` | Whether the current hydrostatic-F inventory contains the analytical barotropic column |

Normalized profiles equal the unscaled common-family profiles, after applying `profileRelation`, divided by the lazily evaluated `normalizationFactors()`, column by column. Wave families already define F=DG-prime; no wrapper conversion is applied. Every field in a column has the same divisor. For zero-APV bases, `boundaryResponseMatrix` composes orthogonalization, normalization and selection from canonical endpoints into current columns. Signed metric signatures remain separate. `basis.metadata` is reserved for caller extensions and is not the authority for these records.

The physical surface convention and generalized-energy coefficients are independent inputs. APV and mean-density modes use the coefficients to define their pairing and endpoint coordinates. The canonical zero-APV response profiles use their active endpoint identities. For neutral DD labels the optional coefficients are provenance only: G vanishes at both ends of every interval, including physical boundaries. Interval derivative values use the inside trace; zero extension does not imply a globally classical derivative across a support endpoint.

`constructionRecord.requested.solver` preserves the backend request. `resolved` identifies the actual backend/source class and available solver settings, including the configured native representation when a source retains it. Analytical barotropic and endpoint-empty bypasses explicitly report skipped construction. Empty numerical selections can discard their native representation; `nativeRepresentationAvailable=false` reports that absence without inventing nodes or reconfiguring the solver.

`constructionRecord.selection` keeps the original ladder request and backend selection diagnostics, along with the current `activeSourceColumns`, `columnLabels` and composed `initialColumnIndices` when the mapping is one-to-one. `sourceEigenpairIndices` identifies the backend eigenpairs where recorded, and remains empty where unavailable; `activeSourceColumns` addresses the extended inventory, which also includes the analytical barotropic column when present. Dense zero-APV combinations use an explicit coordinate map rather than fictitious single source indices. `selectColumns` updates the active identities, spectra and normalization divisors while preserving the original construction request and diagnostics. Requested-empty and subsequently selected-empty bases retain their declared operations.

`constructionRecord.evidence.solveAccuracy` and `referenceConvergence` are explicitly `"unverified"` unless a separate qualification supplies evidence. A correct returned count, a native resolution and an analytical backend name do not substitute for that evidence. Source diagnostics are retained separately.

## Horizontal identity

| `horizontalIdentity.kind` | Families | `k` | `kSquared` |
| --- | --- | --- | --- |
| `"fixedWavenumber"` | Fixed-k waves, inertial auxiliary modes, generalized-enstrophy modes, zero-APV responses and zero-depth waves | One nonnegative solved k per page | Empty |
| `"perColumn"` | Fixed-frequency waves | Empty | One signed horizontal root per column |
| `"independent"` | Hydrostatic F/G, APV, mean-density anomaly and neutral labels | Empty | Empty |
| `"unspecified"` | Developer-supplied scalar bases without a horizontal declaration | Empty | Empty |

The genuine inertial construction has `k=0`. An independent basis does not acquire a fictitious zero wavenumber.

For fixed-frequency modes, the source eigenvalue is the inverse equivalent depth, $$\lambda_j=1/h_j$$. Horizontal dispersion is

$$k_j^2=\frac{(\omega-|f_0|)(\omega+|f_0|)}{g}\lambda_j.$$

Positive roots propagate; negative roots are evanescent. The sign of h alone is insufficient, because the frequency factor changes sign below the inertial frequency. Keeping signed `kSquared` also preserves a zero root without taking a reciprocal of h. The inertial-frequency degeneracy uses the separate inertial factory.

An `IMBasisCollection` accepts caller-supplied `kappa` for explicitly independent bases, so the same vertical structure can serve positive horizontal wavenumbers. Fixed-k pages must match their solved k exactly. Per-column roots and unspecified developer bases reject assignment of a scalar page kappa. With omitted `kappa`, the collection infers a row only when every mapped page has a known scalar wavenumber; otherwise `collection.kappa` remains empty. Each basis's horizontal identity is still available in `collection.metadata`.

```matlab
background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
problem = IMInternalModes.geostrophicAPVModes(background,g0=Inf,gd=Inf);
apv = IMSolverAnalytical().solveEVP(problem,nModes=4);
pages = IMBasisCollection({apv},kappa=[1e-3 2e-3 1e-3],basisIndex=[1 1 1]);
```

## Shared capability authority

`basis.capabilities` is a structure, replacing the former variable-name row in this prerelease API. Its canonical `variables` and parallel `derivativeOrders` cells describe evaluable channels. Evaluation requires a canonical variable and a separate derivative order. Convenience methods such as `Gz(z)` remain available where supported; derivative aliases are not evaluator variables. Supported derivatives are physical vertical derivatives, with interval trace conventions explicit. Scalar evaluation, collection evaluation, chunked evaluation and product factors consume this same record for common bases. Unsupported requests fail before profile evaluation, including on an empty inventory.

```matlab
description = apv.capabilities;
description.variables
description.derivativeOrders
Fz = apv.evaluate([-900;-500;-100],variable="F",derivativeOrder=1);
```

`spectrumNames` distinguishes available spectra from absent ones, including empty spectra. `normalizationNames` lists structurally implemented rules. A retained null or nonfinite column can still make a particular normalization inadmissible; this is reported rather than dropping that column. A scientific pairing can exist before a particular backend implements its normalization operation. Specialized families declare their implemented variables/derivatives explicitly and do not gain F-prime merely by naming an F profile.

`basis.withNormalization(name)` returns an independent basis value using the same continuous representation and selected columns, with an updated active rule and lazily evaluated divisors and signatures. It preserves the original construction requests and evidence. See [changing normalization](BasisNormalization.md) for admissibility, selection invariance and the corresponding changes to coefficients, physical duals and recovery operators.

Zero-APV bases also expose [physical coefficient forms and basis choices](BasisRotations.md). Form getters use current normalized/selected coordinates. Canonical endpoint labels identify source rows; orthogonal column labels identify their linear combinations. Collection metadata preserves both sets and reports that the response map is available on its stored basis; it does not force lazy normalization to copy the map. Its `columnKind="endpoint"` identifies the finite endpoint sector and its fixed-count assessment semantics, rather than asserting that each rotated column belongs to one endpoint.

`projectionOperations` records conditional physical adapters and their observation requirements. Hydrostatic F/G, APV F/G, MDA G and the prescribed neutral APV dual have distinct recipes; see [scalar projection](CommonScalarProjection.md). A defining EVP pairing alone does not promise a full physical state projector.

`sampledRecovery` describes synthesis channels and the requirements for preparing numerical recovery. It does not claim rank, noise control, completeness or recovery success for any observation grid. Fitting evidence, observation/noise factors, activity declarations and measured rank belong to the prepared transform.

## Grid and recovery ownership

The scientific EVP, solver coordinate, observation grid and numerical recovery are separate choices. A provider-designed root grid is one observation geometry; an externally imposed model grid is another. Both use the same continuous identities and can feed the [shared recovery engine](SharedRecovery.md). A diagnostic application may evaluate the required channels on its saved-output grid and pass its assembled observation matrices and factors to `IMDiscreteTransform.fromObservationMatrices`. Boundary traces, derivatives, missing channels and particle labels must be represented by the actual observation model.

This contract does not choose an observation grid, qualify an inverse or define WVM polarization. Model-specific channel assembly and a broader WVM port remain separate work. See [mode creation](ModeCreationAPI.md), [vertical grids](VerticalGrids.md) and [hydrostatic root transforms](HydrostaticRootTransform.md).
