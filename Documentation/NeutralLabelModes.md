# Neutral-label modes in V2

`IMNeutralLabelModesBasis` represents independent stationary displacement labels supported on exact neutral intervals. Their profiles vanish at both ends of each interval, including physical boundaries, and have zero physical and generalized energy. Their amplitudes carry units of length; they have no equivalent depth and are separate from the paired geostrophic F/G families. The implementation follows the manuscript's fixed-endpoint neutral-label vertical structure and APV projection equations. Endpoint particle labels are held fixed; independent boundary displacement anomalies belong to other mode families.

The shared-background creation interface is `IMModes.neutralLabelModes(background,neutralIntervals=...,nModesPerInterval=...)`; it returns an `IMModeBasis` with the same profiles and prescribed APV projector. See the [mode-creation guide](ModeCreationAPI.md) for common identities and normalization names. The direct constructor below remains available.

```matlab
basis = IMNeutralLabelModesBasis([-1000 -850; -700 -400; -250 0],zDomain=[-1000 0],nModes=[2 3 4]);
G = basis.G(z);
q = basis.APV(z);
projection = basis.projectionOnGrid(z,weights);
amplitudes = projection.project(APV);
```

The caller must provide ordered, strictly separated maximal intervals where the reference stratification is exactly neutral. The class neither samples nor infers neutral intervals from density or N2. Merge touching intervals before construction. `nModes` is a nonnegative count of positive-wavenumber modes, scalar for every interval or one count per interval. Zero counts return no columns, including for a fully neutral column. Columns are ordered by interval, then by positive `localModeNumber`; `intervalIndex` and `columnLabels` preserve that identity independently of array indices. The direct basis reports the normalization name `"depthRMSG"`.

With $$D = z_\mathrm{surface} - z_\mathrm{bottom}$$ and $$I = [a,b]$$ of length $$L$$, the dimensionless profiles obey

$$G_{zz} = -m_{\mathrm n}^2 G, \qquad \frac{1}{D}\int_I G_i G_j\,dz = \delta_{ij}.$$

Both interval endpoints have Dirichlet conditions, $$G(a)=G(b)=0$$. The analytical profiles and strictly positive wavenumbers are

$$G_j(z) = \sqrt{2D/L}\sin\!\left(\frac{j\pi(z-a)}{L}\right), \qquad m_{\mathrm n}^j = \frac{j\pi}{L}, \qquad j=1,2,\ldots.$$

This applies to interior intervals, neutral surface and bottom layers, and a fully neutral column. Every row of `boundaryTypes` is `["D","D"]`. The direct constructor's optional `g0` and `gd` arguments are retained as generalized-energy provenance only; they never select label boundary conditions. The common factory likewise accepts optional `energyBoundary` provenance. Changing supported energy coefficients leaves the label profiles unchanged.

The former Neumann endpoint selection, automatically included constant label and `constantLabelPolicy` option have been removed. Existing callers must remove that option and use the positive sine-mode column identities; endpoint anomalies must be represented separately. A fully neutral column with `nModes=3` now has exactly three columns, and `nModes=0` has none. Zero-column raw bases and common views, including `basis.selectColumns([])`, retain empty evaluation and APV projection shapes through `IMBasisCollection`.

`G`, `Gz`, `Gzz` and `APV` return sample-by-mode matrices, with every profile zero outside its interval. `Gz` and `Gzz` use the trace from inside the neutral interval at interval endpoints. `Gzz` is the classical piecewise second derivative; it omits the distributional interface terms of the zero extension. The zero-extended continuous `G` has APV $$q_j = -f_0 G_{j,z}$$, with no interface delta term. `rotationRate` and `latitude` follow the other V2 APIs, with read-only $$f_0 = 2\Omega\sin\phi$$.

`gramMatrix()` returns the exact depth-normalized G Gram, the identity. `gramMatrix(variable="APV")` returns the exact unscaled potential-enstrophy pairing, without a factor of one half:

$$H_{Z,ij} = \int_I q_i^* q_j\,dz = f_0^2 D(m_{\mathrm n}^j)^2\delta_{ij}.$$

Every retained column has positive continuous APV norm when $$f_0\ne 0$$. The factor of one half in physical potential enstrophy gives the norm $$f_0^2D(m_{\mathrm n}^j)^2/2$$. No constant sector needs to be excluded; recovery from particular sampled observations still depends on the supplied rule:

```matlab
basis = IMNeutralLabelModesBasis([-1000 0],zDomain=[-1000 0],nModes=3);
z = -1000+((1:256).'-0.5)*(1000/256);
weights = (1000/256)*ones(size(z));
projection = basis.projectionOnGrid(z,weights);
amplitudes = projection.project(basis.APV(z)*[1;-2;.5]);
```

The projector accepts sampled APV or an APV source directly. For a source state, assemble $$ik\mathcal S_v - i\ell\mathcal S_u - f_0\partial_z\mathcal S_\eta$$ before projection. For an ocean state, provide its labeled APV. The formula is independent of horizontal wavenumber and also applies at zero horizontal wavenumber:

$$A_j = \frac{\sum_r w_r q_j(z_r)^*\mathrm{APV}(z_r)}{f_0^2 D(m_{\mathrm n}^j)^2}.$$

For a displacement source $$s(z)$$, integration by parts gives the equivalent continuous functional

$$\mathcal G_j[s] = \frac{1}{D}\int_a^b G_j(z)s(z)\,dz + \frac{[G'_j(z)s(z)]_a^b}{D(m_{\mathrm n}^j)^2}.$$

The bracket evaluates the source at both interval ends using the neutral-side derivative traces. Endpoint source values need not vanish. Omitting this term changes the physical projector. For any affine source, including the free-surface stretching profile, the volume and endpoint terms cancel. Also $$\int_I q_j\,dz = -f_0[G_j]_a^b=0$$, so constant APV on each neutral interval is orthogonal to every label.


Supply weights for integration in physical `dz`. Interval-local quadrature should resolve each interval and avoid integrating the jump in APV across the stratified side of an interface. For example, concatenate increasing midpoint grids from the separated intervals, using each interval's own midpoint weights. Including stratified points is allowed, but assigning a single trapezoidal weight across an interface can misrepresent the one-sided derivative trace.

The projector divides each pairing by its exact analytical APV norm. It does not fit a sampled Gram inverse, so inaccurate or aliased quadrature remains visible in projected coefficients. It returns an `IMProjection` prescribed dual with a row-scaled identity coefficient system, the physical APV Gram as `majorantGramMatrix`, and the physical Gram recorded in `provenance.continuousAPVGramMatrix`. `projectPhysical` and `project` apply the same physical dual; `recovery` is empty. A caller can separately request sample recovery from `IMDiscreteTransform.fromObservationMatrices`, with an explicit observation/fitting norm. That is a different operation and does not improve the physical quadrature.

For common bases the physical provenance follows the selected columns and their normalization, in the same order as the returned amplitudes. Rank and conditioning describe the identity coefficient system, not grid recoverability; Gram and round-trip assessment capabilities are false. Assess modal recovery explicitly with `projection.forwardMatrix*basis.APV(z)` against the appropriate identity columns. At `f0==0`, APV vanishes for all labels and projection throws `IMNeutralLabelModesBasis:ZeroRotationProjection`; evaluation and exact G normalization remain available.

`IMBasisCollection` accepts this value basis and preserves its label columns. It evaluates `variable="G"` with `derivativeOrder=0,1,2`, and `variable="APV"` at order zero. `collection.projectionOnGrid(z,weights)` defaults to APV for neutral labels and delegates to the same continuous projector; `variable="APV"` and explicit `columns` are supported. Generic Galerkin recipes, generalized-energy projection and automatic collection leakage/product assessment are unavailable for this family. The returned `IMProjection` retains its generic error methods for callers supplying explicit reference coefficients and norms.

Neutral-label amplitudes require independent particle-label information. Velocity and density observations alone do not determine them. The zero-energy statement concerns the fixed linear neutral sector: deleting one of its components from a finite-amplitude displacement can change the reconstructed density when parcel labels cross a plateau boundary. This API does not perform joint state recovery or replace the transport of a chosen label field.
