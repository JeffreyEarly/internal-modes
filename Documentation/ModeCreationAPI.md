# Creating mode families in V2

`IMModes` creates continuous bases from one read-only `IMModeBackground`. Each factory returns an `IMModeBasis` with the same scientific identity, evaluation and normalization vocabulary. Mode requests count the regular ladders; every admissible finite special sector is additional and is returned automatically.

```matlab
background = IMModeBackground.constantStratification(N0=5.2e-3,zDomain=[-1000 0],surfaceCondition="freeSurface");
solver = IMSolverSpectral(nEVP=128,coordinateKind="bandwidth");
k = 2e-4;
waves = IMModes.waveModesAtWavenumber(background,k,solver=solver,nInternalModes=12,nSubinertialModes=0);
z = linspace(-1000,0,201).';
F = waves.F(z);
G = waves.G(z);
waves.modeInfo
```

This free-surface request returns thirteen columns: the surface-gravity mode and twelve internal modes. A rigid-lid background has no surface-gravity mode. `nInternalModes=0,nSubinertialModes=0` requests only the admissible special sector, which can be empty.

## Shared physics and independent solver choices

The background owns `N2`, `zDomain=[bottom surface]`, `rotationRate`, `latitude`, `g` and `surfaceCondition`. The default surface convention is `"freeSurface"`. Its derived quantities are the positive depth `D=diff(zDomain)` and signed Coriolis frequency `f0=2*rotationRate*sind(latitude)`.

For an arbitrary profile, use the constructor directly:

```matlab
background = IMModeBackground(N2=@(z) 1e-4*exp(z/1000),zDomain=[-4000 0],latitude=33,surfaceCondition="rigidLid");
```

`constantStratification` and `exponentialStratification` also record the profile parameters needed by an analytical backend. For an exponential background, `N0` is defined at `zDomain(2)`. The exact exponential backend currently requires that endpoint to be zero. An analytical solver reads the background's recorded source; it owns no independent stratification that could disagree with the physical problem.

```matlab
constant = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
exact = IMModes.waveModesAtWavenumber(constant,k,solver=IMSolverAnalytical(),nInternalModes=12,nSubinertialModes=0);
```

Numerical resolution and coordinate choices belong to the numerical solver. `IMSolverSpectral` accepts `"z"`, classical `"wkb"`, `"bandwidth"` and `"density"` coordinates. Classical WKB uses the monitor $$N$$ and requires strictly positive stratification. Bandwidth uses $$\sqrt{B}$$ with the rotation-aware bandwidth defined by the solver. Neither coordinate changes the physical rotation parameters or adds nontraditional wave dynamics. In particular, `rotationRate=0` changes the physical problem and is not a switch for selecting classical WKB.

The existing finite-difference solver can also be passed to the factories:

```matlab
solver = IMSolverFiniteDifference(z=linspace(-1000,0,257).');
waves = IMModes.waveModesAtWavenumber(constant,k,solver=solver,nInternalModes=8,nSubinertialModes=0);
```

Generalized-energy families additionally require explicit `IMGeneralizedEnergyBoundary(g0=...,gd=...)` coefficients. These describe the basis pairing separately from the physical surface convention. Every finite coefficient, including zero, activates its endpoint coordinate; positive infinity removes it. Negative infinity and NaN are rejected.

## Family catalog and count semantics

| Factory | Requested regular counts | Additional modes returned automatically | Default normalization |
| --- | --- | --- | --- |
| `waveModesAtWavenumber(background,k,...)` | `nInternalModes`, `nSubinertialModes` | Free-surface gravity mode | `"waveEnergy"` |
| `waveModesAtFrequency(background,omega,...)` | `nPropagatingModes`, `nEvanescentModes` | Free-surface mode | `"depthRMSF"` |
| `inertialModes(background,...)` | `nPositiveModes`, `nNegativeModes` of the auxiliary vertical problem | External auxiliary mode for a free surface | `"depthRMSF"` |
| `hydrostaticFModes(background,...)` | `nInternalModes` | External or rigid-lid barotropic sector | `"depthRMSF"` |
| `hydrostaticGModes(background,...)` | `nInternalModes` | Free-surface external mode | `"depthRMSF"` |
| `geostrophicAPVModes(background,...)` | `nPositiveModes` | Every negative and zero eigenvalue mode | `"depthRMSF"` |
| `meanDensityAnomalyModes(background,...)` | `nPositiveModes` | Every negative mode and the admissible constant mode | `"generalizedEnergy"` |
| `geostrophicGeneralizedPotentialEnstrophyModes(background,k,...)` | `nModes` | None | `"generalizedPotentialEnstrophy"` |
| `geostrophicZeroAPVModesAtWavenumber(background,k,...)` | None | Every active surface/bottom response | `"endpointResponse"` |
| `neutralLabelModes(background,...)` | `nModesPerInterval` | Fully neutral constant label when admissible | `"depthRMSG"` |
| `zeroDepthWaveModesAtWavenumber(background,k,...)` | `nModesPerInterval` | None; this explicitly truncates an infinite kernel | `"waveEnergy"` |

All counts are nonnegative integers, and special modes never consume a requested count. Per-interval counts accept a scalar broadcast or one count per declared interval. Numerical requests are exact: insufficient availability throws `IMModes:InsufficientModes` with requested and available counts. Increase solver resolution when the discretization is too small; a physically absent branch cannot be obtained by refining it. A returned count does not certify convergence or recoverability on a separate sampling grid.

Fixed-wavenumber waves are ordered surface first, then internal modes by decreasing frequency, then subinertial modes by increasing frequency. APV and mean-density modes are ordered negative eigenvalues, zero/constant, then the positive ladder. Zero eigenvalue APV modes have infinite equivalent depth; they are distinct from exact zero-depth waves.

`waveModesAtWavenumbers(background,k,...)` returns an `IMBasisCollection`, preserving requested page order and repeated wavenumbers. Counts may be scalar or one value per requested page; repeated wavenumbers must have identical counts and are solved once. Pages may retain different column counts. Zero horizontal wavenumber uses the separate inertial factory.

For example, create the balanced families from the same physical column:

```matlab
background = IMModeBackground.constantStratification(N0=.01,zDomain=[-1000 0]);
solver = IMSolverSpectral(nEVP=128,coordinateKind="z");
energyBoundary = IMGeneralizedEnergyBoundary(g0=-.02,gd=.01);
apv = IMModes.geostrophicAPVModes(background,solver=solver,energyBoundary=energyBoundary,nPositiveModes=12);
mda = IMModes.meanDensityAnomalyModes(background,solver=solver,energyBoundary=energyBoundary,nPositiveModes=8);
zeroAPV = IMModes.geostrophicZeroAPVModesAtWavenumber(background,2e-4,solver=solver,energyBoundary=energyBoundary);
```

The zero-APV basis contains both endpoint responses because both coefficients are finite. Their numerical values select the generalized-energy pairing, while the canonical unit-response profiles depend on the physical background and the active endpoint identities.

At a fixed frequency, the requested propagating and evanescent counts exclude the automatic surface mode. Use `inertialModes` at `omega=abs(background.f0)`; the fixed-frequency horizontal dispersion is degenerate there. Inertial modes all have physical frequency $$|f_0|$$. Their positive and negative eigendepths label an auxiliary vertical basis, and their G profile is auxiliary rather than a displacement mode.

The analytical backend supports only families implemented for the declared profile. Arbitrary profiles have no declared exact backend, and mean-density and generalized-potential-enstrophy families currently require numerical solvers. Unsupported combinations produce an explicit error. Generalized-potential-enstrophy creation currently requires `surfaceCondition="freeSurface"` and accepts its existing `alpha0`/`alphaD` choices.

## Profile convention and normalization

Wave, fixed-frequency wave, inertial auxiliary and exact zero-depth wave profiles use one definition:

$$F_j = D G_j', \qquad D = z_{\mathrm{surface}} - z_{\mathrm{bottom}}.$$

Normalization multiplies every field in a column by the same positive factor and preserves this relation. APV, zero-APV, hydrostatic and mean-density families retain their defining coupled or integral profile relations. Neutral labels have no F profile.

| Name | Meaning |
| --- | --- |
| `"depthRMSF"` | $$D^{-1}\int |F_j|^2\,dz = 1$$ |
| `"depthRMSG"` | $$D^{-1}\int |G_j|^2\,dz = 1$$ |
| `"waveEnergy"` | $$D^{-1}\int [|F_j|^2 + (kD)^2|G_j|^2]\,dz = 1$$ for the fixed-wavenumber wave convention |
| `"generalizedEnergy"` | Unit magnitude under the family's signed generalized-energy G pairing |
| `"generalizedPotentialEnstrophy"` | $$D^{-1}H_Z(F_j,F_j)=1$$ for the family's prescribed generalized-potential-enstrophy pairing $$H_Z$$ |
| `"endpointResponse"` | Canonical unit response at the associated endpoint |
| `"maxAbsF"` | $$\max_z |F_j(z)| = 1$$ |
| `"maxAbsG"` | $$\max_z |G_j(z)| = 1$$ |

The names have the same meaning wherever the required profiles or pairing exist. Unsupported names or family pairings produce an explicit error. Signed normalization preserves its signs in `basis.normSignatures`. A genuinely null generalized-energy mode cannot be assigned unit generalized energy; construction reports the null normalization rather than dropping the mode. Endpoint rotations remain coordinate changes in the existing lower-level zero-APV API.

`waveEnergy` reduces to `depthRMSF` at zero horizontal wavenumber. Fixed-frequency creation defaults to `depthRMSF`, because each returned root has its own horizontal dispersion. The new factory interface does not reinterpret legacy names such as `"kConstant"` or `"geostrophic"` as these different definitions.

## Neutral labels and the constant kernel

Neutral intervals are explicit declarations of exact $$N^2 = 0$$, ordered and strictly separated. They are never inferred from sampled stratification. Merge touching intervals first. Boundary conditions are Dirichlet at joins to stratified water. At a physical endpoint they are Neumann only when its generalized-energy coefficient is zero, and Dirichlet otherwise.

```matlab
background = IMModeBackground(N2=@(z) zeros(size(z)),zDomain=[-1000 0]);
energyBoundary = IMGeneralizedEnergyBoundary(g0=0,gd=0);
labels = IMModes.neutralLabelModes(background,neutralIntervals=[-1000 0],energyBoundary=energyBoundary,nModesPerInterval=3);
```

This returns four columns: the constant label $$G=1$$, followed by three positive-wavenumber cosine modes. The constant has local mode number zero, does not consume `nModesPerInterval`, and has zero APV. `nModesPerInterval=0` returns that constant alone. Neutral labels have neither equivalent depths nor EVP eigenvalues.

APV cannot determine the constant amplitude. A projection request containing it throws `IMNeutralLabelModesBasis:UnobservableConstantLabel`. Select the observable columns explicitly and supply the constant amplitude from independent particle-label information:

```matlab
z = -1000+((1:256).'-0.5)*(1000/256);
weights = (1000/256)*ones(size(z));
columns = find(labels.modeInfo.branch == "neutralLabel").';
projection = labels.projectionOnGrid(z,weights,columns=columns);
sampledAPV = labels.APV(z)*[5;1;-2;.5];
positiveAmplitudes = projection.project(sampledAPV);
```

The projector uses the prescribed continuous APV dual and fixed physical-depth quadrature. Quadrature error remains visible in the coefficients. Zero rotation makes every neutral APV column unobservable and is rejected by projection. See [neutral-label profiles and projection](NeutralLabelModes.md) for the exact formulas and endpoint traces.

## Exact zero-depth waves

Exact $$h=0$$ waves are supported on explicitly declared critical intervals where $$N^2=f_0^2$$. Their kernel can be infinite dimensional, so this family requires its own finite auxiliary-basis request. Small numerical eigendepths are not classified as this exact kernel.

```matlab
f0 = 2*7.2921e-5*sind(33);
background = IMModeBackground(N2=@(z) f0^2*ones(size(z)),zDomain=[-1000 0]);
critical = IMModes.zeroDepthWaveModesAtWavenumber(background,2e-4,criticalIntervals=[-1000 0],nModesPerInterval=4);
```

The implementation uses interval-local Dirichlet sine modes, zero outside their intervals, with $$F=DG'$$ and physical frequency $$|f_0|$$. The caller owns the declaration that the intervals are exact critical regions. This auxiliary truncation is distinct from automatic inclusion of finite special sectors.

## Returned basis and compatibility

`basis.nColumns` reports the full returned inventory. `basis.modeInfo` is a table with `columnLabel`, `branch` and `branchIndex`, plus relevant family metadata such as interval and local mode identity. `columnLabels` preserves identity independently of eigenvalue sign or array position. `F`, `G`, `Gz`, `Gzz` and, when available, `APV` return sample-by-column matrices. Inspect `basis.capabilities` before requesting a channel. `h` and `eigenvalues` report unavailable quantities explicitly for families that do not define them.

The creation layer provides the neutral-label APV projector above. A generic shared wave projection/transform interface is not yet provided by `IMModeBasis`; the existing lower-level basis APIs remain available for their supported transforms and assessments.

`IMInternalModes`, canonical boundary objects, existing numerical/analytical bases and their direct solver entry points remain available. Their original F/G equations, normalization names and count semantics are preserved. In particular, a direct lower-level wave basis retains its existing profile convention; the standard $$F=DG'$$ convention is applied by the new `IMModes`/`IMModeBasis` creation layer. The direct neutral-label constructor now also includes its admissible constant by default; explicit legacy `constantLabelPolicy="reject"` and `"exclude"` remain available there. Existing bulk APIs continue to use their documented lower-level contracts.
