# Common scalar physical projection and recovery

The common continuous basis now supplies scalar physical pairings for hydrostatic F/G, geostrophic APV F/G and mean-density-anomaly (MDA) G. Neutral-label APV retains its prescribed continuous dual. This is the scalar adapter increment of [#48](https://github.com/JeffreyEarly/internal-modes/issues/48). Complete physical state and source projectors remain WVM responsibilities.

## Preparation and application

```matlab
boundary = IMGeneralizedEnergyBoundary(g0=0.2,gd=-0.1);
problem = IMInternalModes.geostrophicAPVModes(background,g0=boundary.g0,gd=boundary.gd);
basis = IMSolverSpectral(nEVP=128).solveEVP(problem,nModes=8);
physical = basis.projectionOnGrid(z,w,variable="G");
aPhysical = physical.projectPhysical(sourceG);

transform = basis.scalarTransform(z,w,variable="G");
aRecovered = transform.transformForward(sourceG);
reconstructedG = transform.transformBack(aRecovered);
```

The point column `z` is caller supplied and contains at least two strictly increasing points. Its aligned physical-depth quadrature column `w` may be signed and must contain a nonzero weight. Preparation preserves both. They may come from another model. An `IMVerticalGrid` is useful for designing and storing rules, but is not required by this scalar adapter. Interior-only quadrature is allowed when the physical pairing has no endpoint terms. When it does, supply the required endpoints as additional scalar observations; their volume weights may be zero.

Specify the scalar `variable` explicitly when preparing a common physical adapter. Supported channels and their required observations are listed below.

Scientific selection and normalization belong to the continuous basis. Select or normalize it before preparing a transform. `projectionOnGrid(...,columns=indices)` constructs a physical projection in those retained coordinates; an explicit `columns=[]` selects no columns. Omitting columns selects all. Collection projection follows the same convention.

## Four separate matrices

For scalar synthesis \(y=Aa\), physical observation metric \(W\), and the continuous signed coefficient form \(H\), the record retains

\[
P=A^*W,\qquad B=H^{-1}P.
\]

`samplePairingMatrix` is P, `targetGramMatrix` is H, and `physicalDualMatrix` is B. `projectPhysical` applies B to independent sampled fields. The positive `majorantGramMatrix` measures coefficient errors; it is formed using absolute physical pairing terms, not by taking the absolute value of H.

The scalar transform separately owns numerical recovery R. Its default is the explicitly named `correctedPhysicalDual` policy:

\[
R=(BA)^{-1}B.
\]

The [shared engine](SharedRecovery.md) prepares that map with scaled rank and overlap checks. B remains unchanged. R corrects the sampled modal overlap; a small recovery round trip does not establish physical quadrature accuracy or accuracy of the continuous reference H. For an independent field outside the retained span, physical coefficients and recovered coefficients generally differ.

To request least squares, declare its fitting factor F or observation factor D:

```matlab
fit = basis.scalarTransform(z,w,variable="G", ...
    policy="sampledLeastSquares",fitFactor=F);
```

F defines the fitting objective. D defines a positive observation/noise norm and supplies F only when F is omitted. Physical quadrature, signed pairings and positive majorants do not silently declare a noise model. The default coefficient scale uses the full active continuous majorant. Recovery policy is independent of grid kind.

A caller can also pass an already prepared `IMDiscreteTransform` using `recovery=prepared`. Its synthesis, labels, activity and physical dual must match exactly. Binding does not prepare or refactor recovery again; physical preparation still constructs its pairing, forms and B. The scalar adapter stores the physical record as `projection` and the generic recovery data on its superclass. Both accept arbitrary trailing field dimensions and complex independent inputs.

Scalar convenience properties mirror the projection record. With DD's prescribed dual, `gramMatrix` and `targetGramMatrix` are its explicitly prescribed coefficient systems (identity here), and `sampledGramRank` is their rank. They do not measure a sampled APV Gram or observation rank. Check `projection.supportsGramAssessment` before interpreting these as physical Gram evidence; use recovery diagnostics for the actual sampled observability. DD's physical APV norm remains its majorant and provenance.

## Family and observation inventory

All expressions below pair the conjugate retained mode with an independent scalar input X. The physical domain is bottom b to surface s. The background supplies N² and gravitational acceleration g.

| Family/channel | Scalar physical pairing or status | Independent source observations |
| --- | --- | --- |
| Hydrostatic F | Interior integral of F* X | F values on the supplied rule |
| Hydrostatic G | Interior integral of (N²/g) G* X; free-surface term Gs* Xs when present | G values, plus surface G for a nonzero surface term |
| APV F | Interior integral of F* X | F values |
| APV G | Interior integral of (N²/g) G* X, plus cs Gs* Xs + (gd/g) Gb* Xb | G values and each endpoint with a nonzero coefficient |
| MDA G | Interior integral of (N²/g) G* X, plus (g0/g) Gs* Xs + (gd/g) Gb* Xb | G values and each endpoint with a nonzero coefficient |
| MDA F | No independent physical scalar projector | Generic sampled recovery may still be prepared from observable F matrices |
| Neutral-label DD APV | Existing exact-norm prescribed continuous APV dual | APV values and quadrature resolving each neutral interval; nonzero f0 |
| Generalized-potential-enstrophy F/G | Scalar physical recipe not bridged here; its G pairing requires independent companion F endpoint traces | Vertical traces are available; the caller must supply the complete physical pairing |
| Wave/inertial/fixed-frequency F/G | Vertical profiles are available; WVM owns complete physical state/source projectors | The auxiliary G EVP pairing alone does not specify physical displacement or complete state/source recovery |
| Geostrophic zero-APV F/G | Coefficient-space energy and boundary forms do not define an independent scalar projector | WVM assembles its paired physical fields and independent endpoint observations |

For APV G, cs is g0/(g+g0) with a free surface and g0/g with a rigid lid. The actual EVP catalog defines the limiting cases: free-surface g0=-g constrains Gs to zero and omits the singular term; free-surface g0=+Inf gives cs=1; rigid-lid g0=+Inf and bottom gd=+Inf omit their terms. Zero coefficients require no endpoint observation. MDA likewise omits inactive infinite endpoint terms.

Endpoint data belong to the independent input. A retained mode's boundary law cannot reconstruct an arbitrary input's missing endpoint value, derivative or companion variable. Value-only recipes reject missing required endpoints and unsupported derivative/companion traces. Adding modes or changing a normalization does not supply those observations.

## Activity, selection and failure

Activity is scientific metadata. The exact F-constant, G-zero sector of rigid-lid hydrostatic modes and literal APV Neumann/Neumann limits is inactive in G, with zero synthesis and coefficient rows; F remains active. A nonzero-norm internal mode that vanishes on the chosen samples remains active and fails numerical recovery if those observations do not resolve it. An active signed-null continuous norm cannot define this physical scalar dual and fails explicitly; it is not converted into an inactive coordinate.

Targets and positive majorants retain their full matrices. Selection and normalization transform those matrices by congruence. Physical preparation solves in the selected coordinates; selecting rows of a previously prepared full-space B generally gives a different map. String `columnLabels` carry common scientific identities. Integer `modeNumber` is optional metadata on the scalar adapter.

Numerical scalar and aligned F/G reference forms integrate products of the represented profiles using oversampled spectral quadrature or per-interval finite-difference reconstruction quadrature. Scalar targets retain off-diagonal entries when the represented columns are not continuously orthogonal. The majorant takes the absolute coefficient of each signed volume or endpoint contribution before summation. These reference forms remain independent of the supplied observation quadrature; a physical dual does not replace its reference Gram with the sampled overlap. Reference-quadrature convergence and eigensolve accuracy remain unverified unless separately established.

Diagnostic recovery failure returns `success=false`, a status and no applicable R. Synthesis and a valid physical B remain available. Calling `transformForward` then throws. `projectionOnGrid` defaults to physical preparation alone. Request `prepareRecovery=true` to expose an additional Galerkin record as `projection.recovery`, then call `recovery.transformForward` or inspect `recovery.forwardMatrix`. `IMProjection` has no public `project` or `forwardMatrix` alias. Use `scalarTransform` to choose another explicit recovery policy. Direct low-level matrix construction with `IMProjection(...)` retains its `prepareRecovery=true` default; pass false there for a physical-only record. A family without a scalar physical pairing can still use the generic matrix recovery interface; the scalar physical adapter does not invent that pairing.

## Reference and qualification boundary

The bridge exposes the existing continuous reference forms and their provenance. Numerical source recipes use solver integration; analytical source recipes use Gauss–Legendre quadrature with `max(128,4*originalColumnCount+32)` points. Selection preserves that original integration budget. Neither is automatically convergence-qualified. Tests compare supported physical forms against independently refined quadrature, including signed endpoint terms and independent inputs, but do not establish a universal retained band or arbitrary-profile convergence.

The aligned `basis.discreteTransform(grid)` remains the rigid-lid hydrostatic F/G adapter. Root-grid construction, external observation-layout qualification, provider grid calculus and WVM-owned physical projector assembly retain their separate scopes. See [verification](CommonScalarProjectionVerification.md) for the evidence and remaining limits of this increment.
