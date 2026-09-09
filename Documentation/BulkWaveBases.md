# Bulk wave bases and bounded evaluation

[Issue #20](https://github.com/JeffreyEarly/internal-modes/issues/20) implements exact bulk computation using the collection contract from [#18](https://github.com/JeffreyEarly/internal-modes/issues/18). InternalModes owns continuous bases and numerical preparation. WVM retains physical polarizations, Fourier interactions, and model-state count decisions.

## Generate a collection

```matlab
N2 = @(z) 1e-4*(1+0.2*exp(z/300));
solver = IMSolverSpectral(nEVP=96);
[collection,costs] = solver.solveWaveModesAtWavenumbers([2e-4 0 1e-4 2e-4],N2=N2,zDomain=[-1000 0],f0=1e-4,nModes=16,nInertialModes=8);
z = linspace(-1000,0,2048).';
F = collection.evaluate(z,variable="F",pages=[4 1 3]);
Gz = collection.evaluate(z,variable="G",derivativeOrder=1,pages=[1 3]);
```

Every distinct requested κ is solved once. `collection.kappa` preserves input order, `basisIndex` maps each request to a stored continuous basis, and every basis retains its actual κ-specific EVP, scientific labels, normalization, and boundary conditions. Positive κ uses the uniform `nModes` count. A zero request requires a separate explicit `nInertialModes`; a solve that cannot return an explicit count errors rather than reducing it. Evaluate heterogeneous wave/inertial counts as separate page selections or with explicitly selected common columns.

The bulk entry supports the built-in spectral and finite-difference solvers. Unknown solver subclasses receive an explicit unsupported-operation error because their preparation may depend on κ. They can continue using independent `solveEVP` calls. Spectral physical, WKB, and density coordinates are covered.

Preparation shares coordinate/grid construction, differentiation matrices, stratification samples for the pencil, and κ-independent assembly. For the wave factory, `A(κ)=A0+κ²*A2` and `B` is shared. Both endpoint rows of `A2` are zero because boundary equations replace those rows. Shared finite-eigenpair processing, scientific selection, basis construction, and orientation preserve scalar behavior. No persistent global cache or approximate neighboring-κ reuse is introduced.

`costs` reports requested/solved κ, mapping, exact reuse counts, preparation count, and setup, shared/page assembly, eigensolve, finalization, collection packaging, and total elapsed time. Shared assembly storage is a component allocation estimate. `solveAccuracy="unverified"` explicitly distinguishes computational provenance from continuous EVP qualification.

## Bound evaluation work

```matlab
F = collection.evaluate(z,variable="F",pages=[4 1 3],sampleChunkSize=256,pageChunkSize=8);
collection.evaluateChunks(z,@consume,variable="F",pages=[4 1 4],sampleChunkSize=256,pageChunkSize=8);

function consume(values,rowIndices,pagePositions,requestedPages)
    % Write or reduce this chunk here; no full result is retained by the API.
    % pagePositions indexes the requested pages vector, so duplicates are unambiguous.
end
```

`evaluate` returns the complete requested `nZ × nColumns × nPages` result and bounds its intermediate work. `evaluateChunks` streams four-argument callbacks and allocates no complete result. Defaults are 1024 sample rows and 16 requested pages per chunk. Selected pages, columns, duplicates, derivatives, and ordering retain the scalar conventions.

Compatible native spectral representations share evaluation/derivative preparation within bounded page groups. Raw native columns are evaluated, EVP-owned diagnostic recovery is applied, and normalization is applied afterward in the same order as scalar calls. Analytical bases, MDA, and custom basis subclasses retain their own evaluators. Numerical and built-in analytical zero-APV evaluators now select source pages before expensive evaluation and preserve rotations and required endpoint identities. Arbitrary custom callbacks still control their own internal allocations; chunk bounds do not claim to control user code.

## Bounded measurements

The construction benchmark uses `nEVP=64`, 8 modes, 16 requested pages and 8 distinct κ, with five measured repetitions after warmup:

| Construction path | Median seconds |
| --- | ---: |
| Independent scalar solve for every requested page | 0.0926 |
| Independent scalar solve for each distinct κ | 0.0680 |
| Bulk solve of the same eight distinct κ | 0.0494 |

This separates exact deduplication from shared preparation: approximately 1.87× versus repeated scalar calls and 1.38× versus already deduplicated scalar calls. In this small case, basis finalization was more expensive than the eigensolve. These observations are not general speedup guarantees. Raw measurements are in [Benchmarks/bulk-construction.json](Benchmarks/bulk-construction.json).

The evaluation benchmark uses 24 distinct/48 requested pages, `nEVP=96`, 16 modes, 2048 rows, and 256-row/8-page chunks. Recorded three-run medians were 0.1704 seconds for scalar F evaluation, 0.0485 seconds for full bulk output, and 0.0522 seconds for streaming. The recorded full outputs matched the scalar result exactly. Full output occupied 12 MiB; the measured largest streaming callback payload was 256 KiB. Raw timings and configuration are in [Benchmarks/bulk-evaluation-final.json](Benchmarks/bulk-evaluation-final.json).

A separate pair of fresh MATLAB processes used 32 distinct/128 requested pages, 4097 rows, 16 modes, and 512-row/8-page chunks. Full evaluation retained 67,125,248 bytes; streaming retained no full result. Maximum resident set size was 803,717,120 versus 744,570,880 bytes; macOS peak footprint was 489,031,432 versus 429,901,528 bytes. These single-process measurements include MATLAB startup and construction and should not be interpreted as precise array allocator peaks. See [Benchmarks/bulk-process-memory.json](Benchmarks/bulk-process-memory.json).

Reproduce construction with `tools/benchmarkBulkWaveModes.m` and evaluation with `tools/benchmarkBulkBasisEvaluation.m`. Keep setup, solve, evaluation, component allocation, and process peak measurements separate.

## Verification and limits

The solver increment passed seven new bulk tests, 43 existing EVP refactor tests, and 19 endpoint cases. The evaluation increment passed six new tests plus 16 numerical and 11 analytical boundary tests. Coverage includes repeated/unordered κ, zero/inertial counts, nondefault/free boundaries, negative/zero labels, three spectral coordinate choices, finite-difference agreement, selected/rotated boundary pages, streaming limits, and custom diagnostic recovery ordering.

The final integrated assessment/collection/bulk batch passed 45 tests. Generated API documentation was refreshed once. The documentation link check, Code Analyzer findings, and package/source-scope checks are recorded in [BulkWaveVerification.md](BulkWaveVerification.md).

This is exact computation and evaluation. It does not assess wave physical residuals, change model-state counts, or introduce a default approximation. The WVM study's documented physical-residual limitation remains relevant. [#19](https://github.com/JeffreyEarly/internal-modes/issues/19) supplies budgeted mixed-product assessment; [#21](https://github.com/JeffreyEarly/internal-modes/issues/21) remains a separate validated-reuse investigation using relevant [#10](https://github.com/JeffreyEarly/internal-modes/issues/10) diagnostics. Provider release/export must precede a WVM dependency update.
