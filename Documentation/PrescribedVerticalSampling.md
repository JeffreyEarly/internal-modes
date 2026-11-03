# Prescribed vertical sampling

`IMVerticalSampling` evaluates an existing continuous `IMModeBasis` at caller-supplied points or over caller-supplied cells. It describes the observations in saved output without choosing a new grid, fitting quadrature weights, changing mode counts or solving another EVP. This is the bounded provider increment of [#60](https://github.com/JeffreyEarly/internal-modes/issues/60).

The sequence is: define the EVP and continuous basis, describe the observations, evaluate their mode matrix, then explicitly prepare numerical recovery. `IMVerticalGrid` remains useful when designing points and quadrature. A sampling layout does not require that object, sorted rows, common F/G locations or endpoint samples. Mode-root points are one possible point layout; they do not change this sequence.

## Observation rows

For the selected profile or classical derivative \(\phi_j\), observation row \(r\) of the synthesis matrix is one of

\[
A_{rj}=\phi_j(z_r),\qquad
A_{rj}=\frac{1}{b_r-a_r}\int_{a_r}^{b_r}\phi_j(z)\,dz,\qquad
A_{rj}=\int_{a_r}^{b_r}\phi_j(z)\,dz.
\]

These are respectively `point`, `cellAverage` and `cellIntegral`. Integrals use physical depth in metres. A cell-center value is not substituted for an average. The basis supplies available variables, aliases and derivative orders; unsupported requests fail.

The immutable `observations` table has one row per observation:

| Column | Meaning |
| --- | --- |
| `observationLabel` | Unique nonempty row identity |
| `variable` | Basis profile or supported alias, such as F, G, APV or Gz |
| `derivativeOrder` | Additional requested derivative order |
| `kind` | `point`, `cellAverage` or `cellIntegral` |
| `zLower`, `zUpper` | Equal for a point; strictly increasing physical cell bounds otherwise |
| `units` | Caller-declared output units; default `unspecified` |

The table constructor supplies `derivativeOrder=0` and `units="unspecified"` when those columns are omitted; the other five columns are required. Additional table columns are preserved as caller metadata. Rows retain their exact ordering. Repeated locations and overlapping cells are allowed when their labels differ. All points and cell bounds must lie inside the basis domain. Cells need a representable interior midpoint; evaluation also rejects subdivisions whose Gauss nodes round onto their endpoints. Units are recorded without conversion or inference. For a cell integral, the caller's declared units must include the depth factor.

## Prepare once, reuse

```matlab
f = IMVerticalSampling.pointValues("F",[-910;-550;-120]);
g = IMVerticalSampling.cellAverages("G",[-970 -700;-620 -260;-200 -30]);
sampling = IMVerticalSampling([f.observations;g.observations]);
[A,record] = sampling.evaluate(basis);
transform = IMDiscreteTransform.fromObservationMatrices(A, ...
    columnLabels=record.columnLabels, ...
    activeColumnMask=true(1,basis.nColumns), ...
    fitFactor=eye(size(A,1)),metadata=record);
a = transform.transformForward(savedObservations);
representedObservations = transform.transformBack(a);
```

`cellIntegrals` is the corresponding convenience constructor for integrated values. Each builder accepts `observationLabels`, `derivativeOrder` and `units`. For more general row mixtures, pass a table with the schema above to the constructor. Use distinct labels when combining blocks with the same variable and observation kind.

The evaluated matrix preserves the basis's column order, normalization and rotations. Its record retains the observation table, scientific and construction provenance, coefficient labels and integration evidence. Pass that record as generic transform metadata. Subsequent forward and backward applications use prepared matrices, including for complex inputs and arbitrary trailing batch dimensions; they do not evaluate modes or integrate cells again.

Activity must be declared for the shared coefficient problem. For example, rigid-lid hydrostatic F/G has a barotropic coefficient whose G profile is zero but whose F profile is constant. It remains active when both channels are observed. G-only observations cannot recover it; the engine reports the deficiency instead of removing that coefficient. No activity mask, fitting objective or observation-noise norm is inferred by the sampling layout.

The example above explicitly chooses an unweighted least-squares objective. It does not declare a noise model. Other fitting factors, whitening factors, coefficient scales and caller-supplied physical duals use the existing [shared recovery contract](SharedRecovery.md). Mixing row units may require an explicit scaled objective. A successful modal round trip establishes finite-span recovery; a field outside that span can agree at every observation while differing elsewhere.

See [the executable staggered-layout example](../ExamplesV2/PrescribedVerticalSampling.m).

## Cell integration evidence

Cell integration evaluates the continuous profile itself. Cells are split at declared neutral or critical interval endpoints, actual finite-difference interpolation knots (including distinct companion-profile knots), and caller-supplied `breakpoints`. Interior Gauss nodes avoid choosing endpoint traces inside an integral. A point at a join follows the basis's existing trace convention.

Classical derivatives retain their basis meaning. They do not acquire distributional interface terms. In particular, integrating a stored finite-difference derivative is not automatically replaced by the difference of endpoint values: the stored derivative and profile may be different interpolated functions. The same caution applies to piecewise neutral-label derivatives at joins.

The default initial Gauss order uses the immutable original column inventory and available solver resolution, with a minimum of 128. Selecting a few high modes does not lower that budget. Successive orders double until every entry meets the requested estimate or `maxOrder` is exhausted. The within-call quadrature rules are reused; there is no persistent cache.

For each cell and column, the estimated error is the sum of absolute differences between successive rules over its pieces. The relative scale is the finer rule's estimate of the integral of the absolute profile, divided by cell width for an average. Thus cancellation in a nearly zero integral does not produce a misleading relative test. The acceptance condition is

\[
e_{rj}\leq\mathrm{absoluteTolerance}_r+\mathrm{relativeTolerance}\,L_{rj}.
\]

Defaults are `relativeTolerance=1e-11`, `absoluteTolerance=0` and `maxOrder=8192`. `absoluteTolerance` may be scalar or one value per observation in the evaluated mode-response units of the synthesis entry, consistent with the supplied basis normalization. These are observation units per coefficient unit; the `units` labels do not convert tolerances. `initialOrder` may be explicitly supplied; `breakpoints` can describe additional known joins. Preparation requires room for at least one doubled rule. Budget exhaustion fails explicitly and returns no sampled matrix. Floating-point range checks also reject nonfinite results, lost positive quadrature measure and complete loss of sampled nonzero profile mass.

The record labels successful integration `estimatedToleranceSatisfied` and retains the per-entry estimated errors, absolute-profile scales, final orders and subdivisions. This is a convergence estimate, not a certified bound for arbitrary evaluator functions. Undeclared narrow features or oscillations may evade paired quadrature rules. Stronger accuracy claims require an independently resolved reference, suitable breakpoints and order controls. Integration evidence also does not qualify the underlying EVP, continuous normalization, physical projection or excluded modal content.

## Responsibility boundary

InternalModes supplies vertical profiles, their point/cell evaluations and generic stable matrix recovery. WVM owns the physical combinations of PV, boundary buoyancy, wave, inertial, mean-density-anomaly and neutral-label coordinates, including horizontal factors, polarizations and source constraints. WVM may assemble a synthesis matrix and optional physical dual from these ingredients and prepare them with the same engine. This increment does not move those physical projector definitions into InternalModes. See [verification](PrescribedVerticalSamplingVerification.md) for the numerical controls and their limits.
