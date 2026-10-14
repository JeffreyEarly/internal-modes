# Neutral-label modes in V2

`IMNeutralLabelModesBasis` represents independent stationary displacement labels supported on exact neutral intervals. They have zero physical and generalized energy under the selected endpoint conditions. Their amplitudes carry units of length; they have no equivalent depth and are separate from the paired geostrophic F/G families. The implementation follows the manuscript's neutral-label vertical structure, generalized-energy endpoint conditions and APV projection equations.

```matlab
basis = IMNeutralLabelModesBasis([-1000 -850; -700 -400; -250 0],zDomain=[-1000 0],nModes=[2 3 4],g0=0,gd=0);
G = basis.G(z);
q = basis.APV(z);
projection = basis.projectionOnGrid(z,weights);
amplitudes = projection.project(APV);
```

The caller must provide ordered, strictly separated maximal intervals where the reference stratification is exactly neutral. The class neither samples nor infers neutral intervals from density or N2. Merge touching intervals before construction. `nModes` is a positive scalar count for every interval or one positive count per interval. Columns are ordered by interval, then by the positive `localModeNumber`; `intervalIndex` and `columnLabels` preserve that identity independently of array indices.

With $$D = z_\mathrm{surface} - z_\mathrm{bottom}$$ and $$I = [a,b]$$ of length $$L$$, the dimensionless profiles obey

$$G_{zz} = -m_{\mathrm n}^2 G, \qquad \frac{1}{D}\int_I G_i G_j\,dz = \delta_{ij}.$$

The lower/upper boundary types and analytical profiles are:

| Boundary types | Profile on the interval | `m_n`, for local mode `j=1,2,...` |
| --- | --- | --- |
| DD | $$\sqrt{2D/L}\sin(m_{\mathrm n}(z-a))$$ | $$j\pi/L$$ |
| DN | $$\sqrt{2D/L}\sin(m_{\mathrm n}(z-a))$$ | $$(j-1/2)\pi/L$$ |
| ND | $$\sqrt{2D/L}\cos(m_{\mathrm n}(z-a))$$ | $$(j-1/2)\pi/L$$ |
| NN | $$\sqrt{2D/L}\cos(m_{\mathrm n}(z-a))$$ | $$j\pi/L$$ |

Every join to stratified water is Dirichlet. A neutral physical surface is Neumann only when `g0==0`; a neutral physical bottom is Neumann only when `gd==0`. All finite nonzero and infinite accelerations select Dirichlet at the corresponding endpoint, enforcing zero generalized-energy coupling. Thus NN occurs only in a fully neutral column with both accelerations zero. That case requires `constantLabelPolicy="exclude"` explicitly: the constant label is fixed to zero because its APV vanishes and the APV projector cannot recover it. Only the positive NN cosine modes are retained.

`G`, `Gz`, `Gzz` and `APV` return sample-by-mode matrices, with every profile zero outside its interval. `Gz` and `Gzz` use the trace from inside the neutral interval at interval endpoints. `Gzz` is the classical piecewise second derivative; it omits the distributional interface terms of the zero extension. The zero-extended continuous `G` has APV $$q_j = -f_0 G_{j,z}$$, with no interface delta term. `rotationRate` and `latitude` follow the other V2 APIs, with read-only $$f_0 = 2\Omega\sin\phi$$.

`gramMatrix()` returns the exact depth-normalized G Gram, the identity. `gramMatrix(variable="APV")` returns the exact unscaled potential-enstrophy pairing, without a factor of one half:

$$H_{Z,ij} = \int_I q_i^* q_j\,dz = f_0^2 D(m_{\mathrm n}^j)^2\delta_{ij}.$$

The projector accepts sampled APV or an APV source directly. For a source state, assemble $$ik\mathcal S_v - i\ell\mathcal S_u - f_0\partial_z\mathcal S_\eta$$ before projection. For an ocean state, provide its labeled APV. The formula is independent of horizontal wavenumber and also applies at zero horizontal wavenumber:

$$A_j = \frac{\sum_r w_r q_j(z_r)^*\mathrm{APV}(z_r)}{f_0^2 D(m_{\mathrm n}^j)^2}.$$

Supply weights for integration in physical `dz`. Interval-local quadrature should resolve each interval and avoid integrating the jump in APV across the stratified side of an interface. For example, concatenate increasing midpoint grids from the separated intervals, using each interval's own midpoint weights. Including stratified points is allowed, but assigning a single trapezoidal weight across an interface can misrepresent the one-sided derivative trace.

The projector divides each pairing by its exact analytical APV norm. It does not fit a sampled Gram inverse, so inaccurate or aliased quadrature remains visible in recovered coefficients. It returns an `IMProjection` prescribed dual with a row-scaled identity coefficient system, the physical APV Gram as `majorantGramMatrix`, and the physical Gram recorded in `provenance.continuousAPVGramMatrix`. Its rank and conditioning describe the identity coefficient system, not grid recoverability; Gram and round-trip assessment capabilities are false. Assess modal recovery explicitly with `projection.forwardMatrix*basis.APV(z)` against the appropriate identity columns. At `f0==0`, APV vanishes for all labels and projection throws `IMNeutralLabelModesBasis:ZeroRotationProjection`; evaluation and exact G normalization remain available.

`IMBasisCollection` accepts this value basis and preserves its label columns. It evaluates `variable="G"` with `derivativeOrder=0,1,2`, and `variable="APV"` at order zero. `collection.projectionOnGrid(z,weights)` defaults to APV for neutral labels and delegates to the same continuous projector; `variable="APV"` and explicit `columns` are supported. Generic Galerkin recipes, generalized-energy projection and automatic collection leakage/product assessment are unavailable for this family. The returned `IMProjection` retains its generic error methods for callers supplying explicit reference coefficients and norms.

Neutral-label amplitudes require independent particle-label information. Velocity and density observations alone do not determine them. The zero-energy statement concerns the fixed linear neutral sector: deleting one of its components from a finite-amplitude displacement can change the reconstructed density when parcel labels cross a plateau boundary. This API does not perform joint state recovery or replace the transport of a chosen label field.
