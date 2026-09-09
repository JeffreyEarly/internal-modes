# Unified basis evaluation and assessment

This is the public contract for [#18](https://github.com/JeffreyEarly/internal-modes/issues/18). InternalModes owns continuous basis generation/evaluation, derivatives, numerical projection, assessment, and solver optimization. WVM supplies physical source channels, polarization factors, valid Fourier interactions, and decisions about model-state counts. No InternalModes API depends on WVM classes.

The API separates four operations:

1. `IMBasisCollection` retains continuous value bases and their scientific identity. Generation supplies this collection; a scalar request is one page.
2. `collection.projectionOnGrid(z,weights,...)` constructs `IMProjection` on exactly the supplied rule and columns. Family-owned pairing descriptors are hidden implementation details, not public recipe objects.
3. `collection.assess(z,weights,...)` or `projection.assess(...)` measures numerical errors without applying tolerances or changing the basis.
4. `checkBasisAssessment(assessment,...)` applies explicit tolerances and returns a separate decision. A rejection or inconclusive reference preserves every requested column.

There is no compatibility constraint on these new interfaces. Existing construction workflows share the numerical projection kernel; their grid design and quadrature fitting remain deliberate operations with separate names. The new fixed-grid entry points do not call fitting or mode selection.

The executable [UnifiedBasisAssessment.m](../ExamplesV2/UnifiedBasisAssessment.m) constructs each family, including scalar, aligned F/G, APV, waves, a separate inertial count, mean-density-anomaly, analytical modes, and both zero-APV boundary implementations.

## Collections and identity

```matlab
collection = IMBasisCollection({basis1,basis2},kappa=[k2 k1 k2],basisIndex=[2 1 2]);
G = collection.evaluate(z,variable="G",pages=[3 1]);
Gz = collection.evaluate(z,variable="G",derivativeOrder=1,pages=2);
```

`bases` stores continuous value objects. `basisIndex` maps each requested page to a stored basis; `sourcePage` selects a page inside existing multi-κ boundary objects. `kappa` preserves requested values/order and must equal any known underlying solved κ. It never substitutes a nearby representative. Metadata preserves the family, normalization snapshot, domain, mode labels or endpoint labels, boundary rotation, variables, and supported derivative orders. Equal κ permits exact sharing only within the same scientific problem, normalization, selection, and discretization context.

Evaluation returns `nZ × nColumns × nSelectedPages`; MATLAB omits trailing singleton dimensions. Pages and explicit column selections preserve caller ordering. Evaluate families with different column counts separately or select an explicit common subset. The initial collection is a correct foundation; shared solver preparation and bounded chunk execution belong to [#20](https://github.com/JeffreyEarly/internal-modes/issues/20).

Numerical mode labels, endpoint identities, κ pages, and frequency signs are distinct. Negative mode numbers are ordinal eigenvalue labels, not permanent boundary-branch identities. Rotated boundary coordinates retain their rotation and endpoint metadata. Boundary coordinates have no automatic retained-prefix policy. Frequency signs remain in caller-supplied interaction identities and do not duplicate a continuous vertical basis. A zero-κ inertial basis/count is independent of the uniform requested wave count.

## Capabilities

| Family or representation | Continuous evaluation | Fixed scalar-variable projection | Assessment limits |
| --- | --- | --- | --- |
| Numerical scalar | `u`, first derivative `uz` | Supplied value-only metric with required endpoints | Gram; positive targets support leakage and products; signed scalar norm policies remain unsupported |
| Aligned F/G, including APV | F/G; first derivative of solved variable | Per-variable signed pairing when all required observations are represented | Gram; positive targets support leakage; signed APV products use a positive majorant for errors |
| Waves and zero-κ inertial modes | F/G; first derivative of solved variable | Projection depends on variable and boundary metric; synthesis does not imply projection | Companion-variable or derivative endpoint observations are explicitly unsupported by value-only projection |
| Mean-density-anomaly modes | F/G; first derivative of solved variable | G projection, diagnostic F synthesis as supported by the family | Coupled quadratic assessment remains explicitly unsupported |
| Analytical aligned modes | Same continuous evaluation contract | Same basis-owned recipe, with no spectral-solver requirement | Analytical functions do not make numerically integrated target Grams exact |
| Numerical/analytical zero-APV boundary bases | F/G and F derivative, with endpoint/page/rotation labels | Automatic scalar F/G projection is explicitly unavailable | Coefficient-space energy/response forms are not scalar sampled projection metrics; no endpoint-prefix selection |

Collection projection and assessment enforce the family capabilities in the table above. Construction checks whether the actual grid includes every required endpoint and whether all required functionals are value-only observations of the selected variable. An unavailable capability raises an actionable error; it is never represented by a successful zero error. Normalization and reference-integration provenance pass into the projection and survive prefix construction.

`IMProjection` is the numerical building block for an explicitly supplied projection recipe, including caller-defined observation spaces. Its basis/metrics are real; projected fields, reference pairings, and coefficient errors may be complex. Physical signed metrics define the sampled Gram solve. The supplied positive majorant defines magnitudes. Gram distortion, round trips, matrix conditioning, and rank have separate meanings; round-trip accuracy alone is not projection or solve qualification.

## Caller-prescribed physical duals

Some source projections use a full-state energy normalization rather than the Gram of a single observed component. WVM wave-source channels are one example. The caller supplies this physical dual explicitly:

```matlab
projection = IMProjection.fromPrescribedDual(samplePairingMatrix,sampleGram,targetGram,majorantGramMatrix=positiveGram,columnLabels=labels,provenance=physicalRecipeProvenance);
assessment = projection.assess(products=preparedProducts,identity=interactionIdentity);
```

`samplePairingMatrix` maps sample values to signed pairings; `sampleGram` is the supplied coefficient-space system used to obtain coefficients, and `targetGram` is the independent continuous system. Both Gram systems and the majorant are real symmetric; pairing matrices and source products may be complex. Provenance is mandatory. Each prefix solves its own supplied system. Required endpoint/source observations must already be represented in the pairing operator and supplied product data.

`projectionKind="prescribedDual"` has no scalar synthesis basis or metric (`sampledBasis` and `metricMatrix` are empty). It reports `sampleCount` and `columnCount` explicitly. Scalar sampled-Gram and round-trip assessment are unsupported, rather than falsely reporting zero error because two supplied normalizations agree. Rank and conditioning of the supplied coefficient system remain available. `projectionKind="galerkin"` retains the basis-derived behavior above. Both use the same product-error and result/decision contract.

Frequency signs can accompany repeated mode numbers in `identity`; distinct coefficient labels must remain unique. `prefixColumnCounts` always counts projection columns. Mapping these to physical wave-mode counts belongs to WVM, including retaining both signs of each selected mode.

## Assessment and decisions

```matlab
assessment = collection.assess(z,weights,page=1,variable="G",columns=1:8,prefixColumnCounts=1:8);
decision = checkBasisAssessment(assessment,gramTolerance=1e-2);
```

Assessments are scalar structs containing `projection`, `identity`, `columnKind`, `prefixColumnCounts`, `measurements`, `coverage`, and `costs`. They carry measurements rather than executing a policy. `measurements` has rows identified by `columnCount` and `quantity` (`gram`, `leakage`, `quadraticAliasing`), with `value`, `status` (`measured`, `inconclusive`, `notRequested`, or `unsupported`), `limitingInputI/J`, `referenceStatus`, `referenceProvenance`, and `examinedCount`. `identity` records the family, variable, page, mapping, normalization, labels, and projection provenance. `coverage` carries inventory metadata; exhaustive and superposition guarantees are always false. `costs` separates projection construction and assessment elapsed seconds. Neither records an independent solve-quality guarantee.

Only requested prefixes are measured. They must increase, stay within the requested band, and include the full requested band. Endpoint-coordinate assessments permit only the complete requested set. Decisions report the unchanged `requestedColumnCount`, raw `requestedBandStatus`, cumulative `status`, and `largestExaminedAcceptedPrefix`. The latter is diagnostic and is not a promise about unexamined prefixes. A previous examined rejection prevents cumulative acceptance; missing or unqualified required measurements are inconclusive. At least one scalar tolerance must be explicitly enabled.

Leakage input is an explicit struct with `sampledValues`, positive `normSquared`, and distinct `columnLabels`. Checks must belong to the same family and page as the output projection; labels identify retained columns to exclude at each prefix. Obtaining additional solved modes is an explicit caller operation.

Product input is an explicit prepared recipe with these fields:

| Field | Meaning |
| --- | --- |
| `sampledValues` | Sample count × product count, after caller-defined physical factors |
| `referencePairings` | Output column count × product count, using the signed continuous pairing |
| `normSquared` | One nonnegative positive-metric product norm square per product |
| `inputColumns` | 2 × product count ordinals in the selected assessment column coordinates (after `columns` reordering), for retained-prefix selection; zero denotes an input outside the output family |
| `inputLabels` | 2 × product count scientific identities; caller can retain page, sign, and interaction labels |
| `referenceStatus` | `qualified`, `inconclusive`, or `unverified`; qualification requires explicit independent evidence |
| `referenceProvenance` | Description of the reference and its quality evidence |

Pairings are supplied instead of full-band reference coefficients because a dense signed target requires a separate continuous Gram solve at each output prefix. Endpoint contributions belong in both the signed pairings and the positive product norms where required. Identically zero products give zero only when their samples and coefficient error are also zero; inconsistent zero norms report infinite error. Norms must never be formed by taking the absolute value of a signed physical quadratic form.

Supplied reference data must be expressed in the assessed basis's normalization and column coordinates. Refining a quadrature rule is separate from independently solving the EVP. A different reference solve must be explicitly matched to the assessed scientific modes; it must never replace them silently. Analytical target Grams currently use a 1024-point trapezoidal integration rule, which is provenance rather than convergence evidence.

[Issue #19](https://github.com/JeffreyEarly/internal-modes/issues/19) extends this prepared-data contract with mixed-basis/page/derivative product inventories, deterministic selection and budgets *before* costly evaluation, and richer family-level coverage and costs. The initial validated selection policy is fixed sparse; targeted additions remain research work.

## Evidence and delivery order

Preserve the frozen cases, signed/trigonometric controls, reference-convergence findings, and measured costs from [WVM #400](https://github.com/JeffreyEarly/wave-vortex-model/issues/400), [study PR #402](https://github.com/JeffreyEarly/wave-vortex-model/pull/402), and [advisory PR #405](https://github.com/JeffreyEarly/wave-vortex-model/pull/405). The documented vertical-momentum/differentiated-pressure residual limitation remains open: sampled-product agreement is not complete model qualification.

[Issue #10](https://github.com/JeffreyEarly/internal-modes/issues/10) owns solve-quality diagnostics and retains its existing milestone commitments. [Issue #21](https://github.com/JeffreyEarly/internal-modes/issues/21) depends on bulk generation and relevant solve-quality diagnostics. It must validate reuse against the requested κ problem and independent references, distinguish accurate-solve acceleration from opt-in approximation, and fall back to a full solve when necessary. No default approximation is authorized.

Release/export the InternalModes provider before WVM updates its dependency to use these capabilities. This API does not change WVM state counts, constructor defaults, or physical source-channel ownership.
