# Creating mode families in V2

The canonical workflow is physical background → family problem → solver → capable family basis. `IMModes` and `IMModeBasis` have been removed. Construction returns the numerical, analytical or specialized family object directly, retaining its applicable diagnostics, normalizations, projections, transforms and boundary-buoyancy basis choices.

```matlab
background = IMBackground.constantStratification(N0=5.2e-3,zDomain=[-1000 0],surfaceCondition="freeSurface");
problem = IMInternalModes.waveModesAtWavenumber(background,k=2e-4);
solver = IMSolverSpectral(nEVP=128,coordinateKind="bandwidth");
waves = solver.solveEVP(problem,nModes=13);
z = linspace(-1000,0,201).';
F = waves.F(z);
G = waves.G(z);
waves.modeInfo
```

Every numerical and analytical `solveEVP` call requires an explicit `nModes`; solver resolution does not choose a scientific mode inventory. Omitting the count raises an error before solving. `nModes=13` requests thirteen total columns. For this free-surface problem that prefix contains one surface-gravity mode and twelve internal modes. Every special sector consumes a requested column. `nModes=0` returns zero columns for either surface convention.

## Shared physics and independent solvers

The read-only background owns `N2`, `rho`, `rho0`, `zDomain=[bottom surface]`, `rotationRate`, `latitude`, `g` and `surfaceCondition`. Rotation and gravity default to Earth rotation, 33 degrees north and 9.81. Every background constructor requires an explicit `surfaceCondition="freeSurface"` or `"rigidLid"`. Its derived quantities are the positive depth `D=diff(zDomain)` and signed Coriolis frequency `f0=2*rotationRate*sind(latitude)`.

Supply exactly one of `N2` or `rho`. Density input sets `rho0` to the density at the upper endpoint of `zDomain`. The returned density and its spectral derivative come from the same Chebfun representation, with $$N^2=-(g/\rho_0)\,\partial_z\rho$$. An explicit `rho0` must agree with the supplied surface density.

Stratification input uses the supplied `N2` evaluator and a Chebfun antiderivative anchored at the surface to derive density. In this case `rho0` defaults to 1025 kg per cubic meter. Both profiles preserve the shape of the queried depths and are independent of solver resolution and coordinates. The input function must describe a fixed profile throughout the background's lifetime.

```matlab
rho = @(z) 1027 - 0.01*z;
background = IMBackground(rho=rho,zDomain=[-1000 0],surfaceCondition="freeSurface");
rho0 = background.rho0;  % 1027 kg/m^3 at z = 0
N2 = background.N2([-1000 0]);
```

```matlab
background = IMBackground(N2=@(z) 1e-4*exp(z/1000),zDomain=[-4000 0],latitude=33,surfaceCondition="rigidLid");
problem = IMInternalModes.geostrophicAPVModes(background,g0=-.02,gd=.01);
solver = IMSolverSpectral(nEVP=128,coordinateKind="z");
apv = solver.solveEVP(problem,nModes=4);
```

Family factories retain their explicit physical-input forms, such as `IMInternalModes.waveModesAtWavenumber(N2=N2,zDomain=[-1000 0],k=k,...)`. With a background, physical inputs and the surface convention must agree with that authoritative column. Without a background, physical factories with a surface choice require an explicit `surfaceBoundary`; their rigid-bottom defaults remain documented. Mean-density-anomaly and generalized-potential-enstrophy factories retain the endpoint laws intrinsic to their named problems rather than accepting an independent surface boundary. The generic mathematical `IMInternalModes` constructor retains its mathematical boundary defaults.

`constantStratification` and `exponentialStratification` backgrounds record the source parameters needed by `IMSolverAnalytical`. For an exponential background, `N0` is defined at `zDomain(2)`; the exact exponential backend currently requires that endpoint to be zero. The analytical solver owns no independent stratification.

```matlab
background = IMBackground.constantStratification(N0=.01,zDomain=[-1000 0],surfaceCondition="freeSurface");
problem = IMInternalModes.waveModesAtWavenumber(background,k=2e-4);
exact = IMSolverAnalytical().solveEVP(problem,nModes=13);
finiteDifference = IMSolverFiniteDifference(z=linspace(-1000,0,257).').solveEVP(problem,nModes=13);
```

An arbitrary profile has no declared exact solution. The profile-independent rigid-lid barotropic sector can nevertheless be constructed exactly, including for a fully neutral column. Mean-density-anomaly and generalized-potential-enstrophy families currently require numerical backends.

Numerical resolution and coordinates belong to the solver. Spectral `"wkb"` means the classical N coordinate and requires strictly positive stratification. `"bandwidth"` uses the rotation-aware squared-frequency bandwidth. `"z"` and `"density"` are also available. Coordinates do not change rotation parameters or add nontraditional wave equations. Use `rotationRate=0` when the physical problem is nonrotating.

## Exact counts and explicit branches

Without `branch`, `nModes` selects an exact prefix of the family's physical ordering. Requests larger than the available inventory raise `IMSolver:InsufficientModes`. A returned count does not certify numerical convergence or recovery on a separate sampling grid.

For explicit inventories, provide one count per requested branch:

```matlab
problem = IMInternalModes.waveModesAtWavenumber(background,k=2e-4);
internal = solver.solveEVP(problem,branch="internal",nModes=3);
mixed = solver.solveEVP(problem,branch=["surface","internal","subinertial"],nModes=[1 3 2]);
```

The first request returns exactly three internal columns and no surface column. The second returns six columns in the requested branch order, provided the physical profile and resolution support those branches. A constant strongly stratified column has no subinertial branch; refining it cannot create one. There is no scalar broadcast across multiple explicit branches. Unknown branch names raise an error even when their requested count is zero.

| Problem factory | Principal inventory or branches | Default normalization |
| --- | --- | --- |
| `waveModesAtWavenumber(background,k=k)` | `surface`, `internal`, `subinertial`; includes auxiliary inertial modes at k=0 | `"waveEnergy"` |
| `waveModesAtFrequency(background,omega=omega)` | `surface`, `propagating`, `evanescent` | `"rmsF"` |
| `hydrostaticFModes(background)` | Rigid-lid `barotropic` plus positive modes, or free-surface external and positive modes | `"rmsF"` |
| `hydrostaticGModes(background)` | Free-surface external and positive modes, or rigid-lid positive modes | `"rmsF"` |
| `geostrophicAPVModes(background,g0=g0,gd=gd)` | `negative`, `zeroEigenvalue`, `positive` | `"rmsF"` |
| `meanDensityAnomalyModes(background,g0=g0,gd=gd)` | `negative`, `constant`, `positive` | `"meanDensityAnomalyEnergy"` |
| `geostrophicGeneralizedPotentialEnstrophyModes(background,k=k)` | Positive generalized-enstrophy ladder | `"generalizedPotentialEnstrophy"` |

Fixed-wavenumber wave prefixes contain the surface mode first, then internal modes by decreasing frequency, then subinertial modes by increasing frequency. APV and mean-density prefixes contain negative eigenvalues, then the zero/constant sector, then positive eigenvalues. A zero APV eigenvalue means `h=Inf`; it is distinct from an exact zero-depth wave.

At `omega=abs(background.f0)`, use `waveModesAtWavenumber(background,k=0)` because fixed-frequency horizontal dispersion is degenerate. All inertial auxiliary columns have physical frequency $$|f_0|$$; their positive and negative equivalent depths label an auxiliary vertical problem, and G is auxiliary rather than a displacement profile.

The bulk API preserves requested page order, repeated wavenumbers and unequal counts. Its required `nModes` is a scalar applied to every page, including zero wavenumber, or a vector aligned with the full `kappa` input. Zero counts are valid; an empty `kappa` is rejected. Each distinct wavenumber is solved once for its largest requested count, with smaller requests represented by exact prefix selections. Identical `(kappa,count)` requests share one stored basis. Each prefix records its own requested count and shared-solve provenance while retaining the original normalization integration budget.

```matlab
collection = solver.solveWaveModesAtWavenumbers(background,[2e-4 0 1e-4 2e-4],nModes=[5 3 7 2]);
F = collection.evaluate(z,variable="F",pages=[4 1],columns=1:2);
```

Numerical bulk construction reuses the prepared pencil; analytical bulk construction reuses exact solutions only for identical wavenumbers. Neither approximates distinct pages.

## Consistent evaluation

Scalar bases and collections use the same explicit variable and derivative options:

```matlab
G = waves.evaluate(z,variable="G");
Gz = waves.evaluate(z,variable="G",derivativeOrder=1,normalization="waveEnergy");
collection = IMBasisCollection({waves});
GzPages = collection.evaluate(z,variable="G",derivativeOrder=1,normalization="waveEnergy");
```

`variable` is required. `derivativeOrder` defaults to zero; derivative names such as `"Gz"` are not variable aliases. Convenient methods such as `F(z)`, `G(z)`, `Gz(z)` and `Gzz(z)` remain available where the family supports them. Collections return `nZ × nColumns × nSelectedPages` values, with selected page and column order preserved. Omitting `pages` selects all pages; explicit `pages=[]` selects none and streams no callbacks. An empty page selection returns `nZ × 0 × 0` when columns are omitted, or `nZ × numel(columns) × 0` with explicit columns. Nonempty page selections with heterogeneous column counts require explicit common columns or separate homogeneous groups.

Ordinary and streamed collection evaluation both accept `normalization`. Omitting it uses each stored basis's default. An override scales the requested profiles and derivatives without changing stored defaults, mode identities or the eigensolve. Collection metadata records each stored default and is unchanged by a per-call override. Unsupported variables, derivative orders and normalization overrides are checked across the selected bases before any streaming callback. Endpoint bases apply the same column scaling to their response maps and physical forms. Orthogonalize the complete endpoint family before selecting columns; normalization remains a separate choice. See [bounded collection evaluation](BulkWaveBases.md).

## Profile convention and lazy normalization

Every public wave construction path, including fixed-frequency, inertial and exact zero-depth waves, uses the manuscript definition

$$F_j = D G_j', \qquad D = z_{\mathrm{surface}} - z_{\mathrm{bottom}}.$$

APV, zero-APV, hydrostatic and mean-density families retain their distinct defining relations. Neutral labels have no F profile. Normalization divides every aligned field in a column by the same positive factor and preserves its diagnostic relations.

```matlab
Fdepth = waves.F(z,normalization="rmsF");
Gdepth = waves.G(z,normalization="rmsF");
waves.normalization = "maxAbsF";
Fmax = waves.F(z);
waves = waves.addNormalization("custom",@(~,~) 2);
```

Custom rules have signature `factor=rule(basis,j)`. Profile factors are computed when requested; numerical and analytical family bases reuse expensive factors without repeating the eigensolve. Reading a spectrum does not evaluate profiles merely to choose a normalization. `selectColumns` retains aligned scientific identities, profiles and family data without solving again.

All built-in choices keep F and G dimensionless. RMS uses the entire physical depth D, including zero extensions outside declared intervals. It refers to the named profile, not a universal velocity or displacement interpretation.

| Name | Meaning where supported |
| --- | --- |
| `"rmsF"`, `"rmsG"` | $$D^{-1}\int |F_j|^2\,dz = 1$$ or $$D^{-1}\int |G_j|^2\,dz = 1$$ |
| `"rms"` | Generic scalar basis: $$D^{-1}\int |u_j|^2\,dz = 1$$ |
| `"maxAbsF"`, `"maxAbsG"` | Unit maximum absolute profile amplitude |
| `"surfaceF"` | $$|F_j(z_\mathrm{surface})|=1$$ |
| `"waveEnergy"` | $$D^{-1}\int [|F_j|^2 + (kD)^2|G_j|^2]\,dz = 1$$ at fixed wavenumber, including k=0 |
| `"generalizedPotentialEnstrophy"` | The complete boundary-augmented F pairing divided by D equals one |
| `"meanDensityAnomalyEnergy"` | The magnitude of $$(\int N^2|G_j|^2\,dz+g_0|G_j(s)|^2+g_d|G_j(b)|^2)/g$$ equals one; the signed state energy is $$g\epsilon_j/2$$ |
| `"endpointResponse"` | Canonical unit endpoint response, or unit Euclidean response-vector length after generalized-energy orthogonalization |

`normalizationNames()` reports supported rules. Neutral labels default to `rmsG`, exact zero-depth waves to `waveEnergy`, endpoint modes to `endpointResponse`, and generic scalar modes to `rms`. Infinite inactive MDA endpoint terms are omitted before arithmetic.

All canonical names are reserved. Custom rules must use a distinct name; replacing that custom name remains supported. V2 rejects legacy aliases and `Normalization` enum values. It offers neither `innerProduct*` variations nor a `generalizedEnergy` normalization. Physical energy matrices and generalized-energy orthogonalization remain available independently of profile scale. Fixed-frequency modes do not offer one fixed-wavenumber `waveEnergy` normalization. A null profile cannot receive its corresponding unit-amplitude normalization and fails explicitly without removing the column.

The returned family objects retain applicable `gramMatrix`, `majorantGramMatrix`, `endpointGramTerms`, derivatives, spectra and transform methods. Neutral labels retain their prescribed APV dual. Endpoint-response families retain their canonical identity, current boundary-response map and quadratic matrices through normalization and selection. A capability present on the returned family is not hidden behind a generic wrapper.

## Explicit zero-APV endpoints

Generalized-energy coefficients define APV problems, pairings and endpoint-basis choices. They do not select canonical endpoint responses.

```matlab
problem = IMGeostrophicZeroAPVModes.atWavenumber(background,2e-4,endpoints="surface");
surface = solver.solveGeostrophicZeroAPVModes(problem);

problem = IMGeostrophicZeroAPVModes.atWavenumber(background,2e-4,endpoints=["surface","bottom"]);
boundary = solver.solveGeostrophicZeroAPVModes(problem);
depth = boundary.orthogonalizeGeneralizedEnergy(g0=-.02,gd=.01);
hbb = depth.boundaryDepth;
M = depth.boundaryResponseMatrix;
Hg = boundary.generalizedEnergyMatrix(g0=-.02,gd=.01);
```

Canonical columns have unit response at the requested endpoint and zero response at the other endpoint. Both endpoint identities are returned in surface/bottom order. `endpoints=strings(1,0)` explicitly requests an empty coordinate inventory. These bases have no ordinary EVP eigenvalues or equivalent-depth property; the generalized-energy orthogonal basis separately exposes signed `boundaryDepth`. Numerical and analytical endpoint families preserve pagewise response, energy and derivative operations. Orthogonalization requires the complete originally requested endpoint family; choose its directions before reducing columns. Negative, zero and repeated boundary depths remain valid. See [boundary-buoyancy basis choices](BasisRotations.md).

## Neutral labels and exact critical intervals

Neutral intervals are caller declarations of exact $$N^2 = 0$$, never inferred from sampled stratification. Declare ordered, strictly separated intervals and merge touching intervals first. Profiles vanish at both interval ends, including physical endpoints. Independent boundary anomalies belong to endpoint-response families.

```matlab
neutral = IMBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0],surfaceCondition="freeSurface");
labels = IMNeutralLabelModesBasis(neutral,neutralIntervals=[-1000 0],nModesPerInterval=3);
z = -1000+((1:256).'-0.5)*(1000/256);
weights = (1000/256)*ones(size(z));
projection = labels.projectionOnGrid(z,weights,variable="APV");
amplitudes = projection.projectPhysical(labels.APV(z)*[1;-2;.5]);
```

This returns exactly three positive-wavenumber DD sine modes. There is no constant label. Specialized interval families require `nModesPerInterval` as an explicit count per interval, scalar broadcast or one count per declared interval; zero retains no columns on that interval. Optional `g0` and `gd` on labels record provenance and do not alter the profiles. Labels have neither ordinary eigenvalues nor equivalent depths. APV projection uses the prescribed continuous dual, including the active normalization, and fixed physical-depth quadrature. Quadrature error remains visible; zero rotation makes APV unobservable and is rejected. See [neutral-label profiles and projection](NeutralLabelModes.md).

Exact $$h=0$$ waves use a separate family on declared critical intervals where $$N^2=f_0^2$$. Small numerical depths are not classified as this exact kernel.

```matlab
f0 = 2*7.2921e-5*sind(33);
critical = IMBackground(N2=@(z) f0^2*ones(size(z)),zDomain=[-1000 0],surfaceCondition="freeSurface");
zeroDepth = IMZeroDepthWaveModesBasis(critical,2e-4,criticalIntervals=[-1000 0],nModesPerInterval=4);
```

The DD sine profiles are zero outside their intervals, satisfy $$F=DG'$$ with one-sided interval derivative traces, and have physical frequency $$|f_0|$$. The caller explicitly truncates the infinite-dimensional kernel. `h` is exactly zero; this family has no finite inverse-depth EVP eigenvalue. Exact analytical profile norms, derivatives, column selection and Gram matrices remain available.

## Migration from the removed creation layer

The earlier facade counted regular modes and appended special sectors. Replace it with an explicit total or branch request:

```matlab
% Before: two negative APV columns were added to two positive columns.
% IMModes.geostrophicAPVModes(background,solver=solver,energyBoundary=energy,nPositiveModes=2)

% After: request four total columns without predicting the negative count.
problem = IMInternalModes.geostrophicAPVModes(background,g0=energy.g0,gd=energy.gd);
apv = solver.solveEVP(problem,nModes=4);
```

Replace `IMModes.geostrophicZeroAPVModesAtWavenumber` with an explicit `IMGeostrophicZeroAPVModes.atWavenumber(...,endpoints=...)` problem and `solveGeostrophicZeroAPVModes`. Construct labels and exact zero-depth waves with their specialized classes. `IMInternalModesDiscreteTransform.hydrostaticRoot` returns a capable native hydrostatic basis; its fixed-grid transform and qualification contracts remain available. See [hydrostatic root transforms](HydrostaticRootTransform.md).

When migrating direct lower-level wave users, update modal coefficients and reconstruction alongside profiles. The provider wave convention is $$F=DG'$$ and its default is `waveEnergy`. WVM preserves its existing model coordinates for finite positive h using $$F_\mathrm{WVM}=\sqrt{h/D}F,\qquad G_\mathrm{WVM}=\sqrt{D/h}G,$$ including the same factors on spectral derivatives and at k=0. WVM owns that conversion; saved model arrays are not rescaled on restore. Mean-density-anomaly scales are unchanged under their new descriptive name. These conversions are not a definition for nonpositive-depth wave branches; each consumer must specify its physical polarization and supported inventory.
