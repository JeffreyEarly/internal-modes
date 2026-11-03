# Shared prepared recovery

V2 separates the scientific basis, physical projection and numerical recovery. This prerelease API intentionally replaces the old scalar `IMDiscreteTransform` constructor. It implements [#57](https://github.com/JeffreyEarly/internal-modes/issues/57); cross-grid qualification and WVM adoption remain separate work. InternalModes supplies vertical ingredients and generic numerical preparation; WVM owns its complete physical state and source projectors.

## Ownership

| Object | Responsibility |
| --- | --- |
| Capable scientific family bases | Continuous scientific families, labels, normalization, evaluation, derivatives and physical recipes |
| `IMVerticalGrid` | Physical points, quadrature weights and grid-design provenance |
| `IMProjection` | Physical pairings, continuous coefficient normalization, sampled/target Gram diagnostics and physical-dual application |
| `IMDiscreteTransform` | Prepared synthesis/recovery matrices, coefficient identity, observation metadata and numerical recovery evidence |
| `IMScalarDiscreteTransform` | Scalar grid/quadrature adapter, with a genuine physical `projection` record |
| `IMInternalModesDiscreteTransform` | Aligned F/G adapter, with separate channel `recovery` and `physicalProjection` records |
| `prepareIMRecovery` | Stateless numerical preparation shared by the adapters |

Creating a continuous basis does not select a grid. Selecting a grid does not select a recovery policy. Preparing from fixed matrices performs no EVP solve, weight fitting or retained-count search. Factory methods may perform these scientific preparation steps explicitly before calling the shared engine.

## Matrix contract

For observations `y` and scientific amplitudes `a`, synthesis is `y=A*a` and recovery is `a=R*y`. `A` is a finite real or complex `sampleCount × columnCount` matrix. Unique string `columnLabels` and an explicit `activeColumnMask` preserve coefficient identity and alignment. Inactive columns must be exactly zero; their recovered rows are exactly zero. An active column that vanishes at the selected observations is unobservable and fails recovery.

Observation rows may stack variables on different points, endpoint coordinates, traces, cell averages or other caller-assembled linear observations. Their count is independent of `length(z)`. The generic object requires no `IMVerticalGrid`; a diagnostic application can use another model's prescribed observations directly. WVM ingestion and physical field combinations remain WVM responsibilities.

The engine keeps four inputs distinct:

- `physicalDual=B`: optional continuous-normalized physical coefficient dual, with aligned coefficient rows. This remains the original physical dual after recovery correction.
- `coefficientScale=C`: nonsingular amplitude map on the declared active coefficient space. A vector supplies a diagonal map; a matrix supplies a full map. Mixing active and inactive sectors is unsupported by this activity contract.
- `observationFactor=D`: nonsingular whitening factor defining the positive observation/noise norm. The scaled synthesis is `As=D*A*C` on active columns.
- `fitFactor=F`: optional rectangular least-squares objective, defaulting to D when D is declared. A nonnegative physical metric can supply a square-root factor including zero weights. F need not equal D.

Neither signed physical forms nor their absolute values define a positive observation norm. Zero volume weights also do not remove endpoint observations. With no D declaration, algebraic recovery can use F or a supplied weak dual, but observation rank/geometry and noise gain remain unavailable. This supports physical adapters without inventing a scientific noise norm. Callers requiring noise or geometric evidence must supply D explicitly.

## Named policies

### Sampled least squares

`policy="sampledLeastSquares"` factors `F*A*C` directly, before any Gram formation. On a full-column-rank fit it constructs the exact left inverse `R=C*(F*A*C)^+*F`. The pseudoinverse notation describes the full-rank least-squares formula; implementation never silently truncates active directions.

When F=D, recovery minimizes observation-noise amplification among exact left inverses in the declared coordinates. The diagnostics report this capability. When F differs from D, they report the actual declared observation-noise gain without claiming optimality. For example, `A=[1;1]`, D=I and `F=[1 0]` gives `[1 0]`; observation-optimal fitting gives `[.5 .5]`. An observation rank of two with a fit rank of one fails the requested fitting policy.

### Corrected physical dual

`policy="correctedPhysicalDual"` preserves the weak equations while correcting sample overlap. With `Bs=C^-1*B*D^-1`, the mathematical equations are `(Bs*As)*L=Bs`, followed by `R=C*L*D`.

The engine equilibrates test equations and uses an orthonormal test-row frame before deciding overlap rank and solving. Invertible left preconditioning changes the representation of the equations, not their recovery map. It separately records observation, test and overlap evidence. Signed pairings can genuinely have deficient overlap despite full sample rank; that fails explicitly. Scaling cannot supply missing information.

### Extended span

`policy="extendedSpan"` prepares an explicitly supplied larger synthesis family using `extendedPolicy`, then keeps `retainedColumnIndices`. Supply aligned extended labels, activity, scales and, for weak recovery, the extended physical dual. Retained A must agree with the selected extended columns. Recovery rejects the represented omitted columns within the measured residual; there is no claim about arbitrary higher modes.

`selectOutputs(columns)` creates a cached view of an already prepared map without refitting. Repeated selections compose to the original extended identities. The original specification and diagnostics remain labeled by their construction span. Selected unscaled round-trip/rejection errors are measured; selected coefficient-norm diagnostics remain unavailable until an explicit retained scale is prepared. Caller metadata is preserved unchanged.

The classical root factory is a specialization: its `truncatedCollocation` channel prepares the full root interpolation span once, then keeps qualified outputs. The physical adapter name and the engine's `extendedSpan` policy describe different levels of the same construction.

## Generic usage and failure

```matlab
A = [1 0; 0 1; 1 1];
transform = IMDiscreteTransform.fromObservationMatrices(A, ...
    columnLabels=["mode:1" "mode:2"], ...
    coefficientScale=[1 1], observationFactor=eye(3), ...
    metadata=struct("layout","three externally prescribed observations"));
a = transform.transformForward(values);
reconstructed = transform.transformBack(a);
```

Both operations preserve arbitrary trailing batch dimensions and use cached matrix multiplication. They do not refactor or repeat scientific preparation.

Default `failureMode="error"` throws for an unresolved requested policy. `failureMode="diagnostic"` returns `success=false`, an explicit `status` and **no forward matrix**. Calling `transformForward` then throws `IMDiscreteTransform:UnavailableRecovery`. Synthesis remains available. There is no automatic tolerance relaxation, policy substitution, pseudoinverse fallback or mode dropping. Empty/inactive sectors retain their declared shapes.

`recoverySpecification` records the resolved factors, coefficient scaling, rank tolerance and extended-span choices independently of caller metadata. Diagnostics distinguish fit rank, observation rank, test rank, overlap conditioning, retained round trip, defining-equation residual, noise gain and extended rejection. Every capability has an availability state. An absent observation norm never produces a successful zero noise measurement. Extremely scaled spectra also retain their scale representation when a dimensional singular value cannot be represented directly.

Finite inputs alone do not guarantee that every intermediate coordinate change or recovery matrix fits floating-point range. Such a preparation fails explicitly; diagnostic mode retains the supplied dual/specification with no R. Full matrix whitening is supported, but no general structured or implicit-whitening memory guarantee is made. Algebraic weak recovery without D avoids allocating an observation-sized identity or projector.

## Physical projection and migration

For a Galerkin physical recipe, `IMProjection` retains raw pairings and the continuous-normalized dual B. `projectPhysical(values)` applies B directly. An explicitly prepared `recovery.transformForward(values)` applies R. These generally differ off the represented modal span. Physical/source projection remains meaningful even if sample recovery is unavailable, provided the physical coefficient normalization is valid. The ambiguous `project` and `forwardMatrix` aliases have been removed from `IMProjection`.

Positive physical metrics supply F to the engine. Signed metrics use corrected weak recovery. A sampled Gram that rounds to deficient rank does not veto a resolvable positive fit: for `A=[1 1;0 1e-8]`, the engine recovers both columns even though `A'*A` loses a direction. The actual physical Gram error and rank remain separate evidence; this does not qualify the quadrature. Scalar/F/G count policies still enforce their physical qualification criteria.

`projectionOnGrid` defaults to `prepareRecovery=false`, creating only the physical projection record. Direct matrix construction with `IMProjection(...)` retains its `prepareRecovery=true` default for low-level assessment callers; pass false there to request only physical preparation. Collection assessment explicitly prepares its required recovery. A prescribed source dual without synthesis remains a physical projection and does not manufacture a recovery object. A singular prescribed coefficient system retains raw pairings with an unavailable physical dual, rather than a least-squares substitute.

An unrepresentable physical coefficient dual similarly has status `unrepresentableCoefficientSystem` and no applicable matrix. This does not veto independent positive-fit sample recovery. Corrected physical recovery requires a representable normalized dual; otherwise construction reports an actionable error, and `prepareRecovery=false` remains available for inspecting the raw physical record.

Migration from the earlier prerelease API:

- Use the common basis `scalarTransform(z,w,variable=...)` factory for supported physical scalar channels, or bind prepared records with `IMScalarDiscreteTransform(physical,recovery,z=z,weights=w,normalization=name)`. String labels and full continuous forms are supported. Legacy basis `discreteTransform` factories continue to prepare the scalar adapter. Scalar physical metrics and Gram evidence live on its `projection` record. See the [family and observation inventory](CommonScalarProjection.md).
- Use `IMDiscreteTransform.fromObservationMatrices` for arbitrary or coupled observations. Generic transforms expose no invented Gram/metric/Parseval fields.
- Use aligned transform `recovery(variable=...)` and `physicalProjection(variable=...)` to inspect numerical and physical records separately. F0/G0 activity remains explicit.
- Replace arbitrary forward-matrix injection with a prepared shared-engine result in provider adapters. Physical names such as `galerkin` and `truncatedCollocation` remain available through `forwardOperator`.

## Qualification boundary

The portable kernel and adapter regressions cover complex algebra, coordinate/unit changes, fit/noise separation, active/inactive sectors and explicit deficiencies. [Verification and cost evidence](SharedRecoveryVerification.md) records the tested revision, independent audits and known baseline failure.

This implementation does not raise retained bands or choose a universal WVM policy. Provider vertical operations belong to #48, provider grid calculus to #50, cross-grid/family qualification to #58 and prescribed-grid fixtures to #60. Complete physical state/source projectors and production WVM integration belong to WVM, tracked downstream in #59. The [neutral-label API](NeutralLabelModes.md) documents the manuscript DD convention separately; local coupled-grid DD research is not presented here as released qualification evidence.
