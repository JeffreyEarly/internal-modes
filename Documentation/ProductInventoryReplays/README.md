# Additional frozen WVM replay verification

Both bounded cases passed the existing `2e-7` absolute comparison gate without changing it. Every signed product/retained-prefix error was compared after conversion to the frozen single precision storage format; every structural-zero flag and nonzero count matched. Prefix maxima additionally compare the new double precision values to the frozen CSV values.

| Frozen fixed case | Products | Nonzero | Zero | Product-prefix comparisons | Max product difference | Max prefix difference |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| cal-exponential-17 | 79,968 | 48,592 | 31,376 | 639,744 | 7.276e-12 | 2.84e-09 |
| withheld-pycnocline-33 | 236,192 | 142,528 | 93,664 | 3,779,072 | 3.638e-12 | 6.591e-09 |

Each verification JSON includes all per-row differences/counts, all prefix values/differences, measured reference discrepancies, original frozen reference diagnostics, stage costs, and runtime class paths. The manifest hashes the executed adapter, scientific and new numerical sources, frozen inputs, complete MAT fixtures, JSON evidence, and execution log. Large fixture MAT files remain temporary and need not be committed.

The exponential case used EVP orders 128/192 and reference quadrature orders 257/513. The withheld pycnocline case used the documented refinements 192/256 and 513/1025. Neither assessed scientific basis was silently replaced.

These records describe the original API before consolidation. Current callers use `IMProductAssessmentPlan` with factor/output structures, and `checkProductAssessment` and `compareProductAssessmentCoverage` on returned assessment structures. Historical source names and hashes below remain unchanged for reproducibility.

## Costs and retained reference limitations

Exponential source preparation took 3.2947 s, adapter preparation including inventory/planning 1.5990 s, and the complete generic assessment call 4.7422 s. Pycnocline took 3.6317 s, 2.5088 s, and 16.5512 s respectively. These are one sequential local replay per case, not a controlled performance comparison; inner engine timing is separately retained in JSON. The adapter constructs the full selected-source cache before handing it to the new API, so these timings do not demonstrate deferred physical source preparation. Cached factor array payloads were 184,650,278 and 855,721,984 bytes; these exclude MATLAB process/allocator overhead.

Maximum measured independent-reference product discrepancies were 5.0876e-6 and 2.7284e-5. The frozen refined studies report stable required references, while their broader all-mode derivative differences remain 4.4031e-5 and 4.9320e-5. The withheld case also retains rejected fixed-family Gram quality. The replay preserves these observations; agreement with stored values does not qualify all derivatives, the eigenproblem, physical residuals, superpositions, or a complete model. Measurement and acceptance policy remain separate.

The scientific solver was loaded from `wvm400-oceankit/InternalModes-2.0.0-beta.1`, and the new `IMProductInventory`/`IMProductPlan`/`IMProductAssessment`/`IMProjection` numerical classes from `internal-modes-product-inventories`. Running from the provider worktree would shadow the pinned scientific solver and is not a valid replay. Use a fresh MATLAB process with `/private/tmp` as the current directory; see the manifest command. Nonfatal startup warnings about the absent personal `Documents/MATLAB` path and WVM package search-path ordering are retained in the log.

Related work: [InternalModes #19](https://github.com/JeffreyEarly/internal-modes/issues/19), [WVM study #400](https://github.com/JeffreyEarly/wave-vortex-model/issues/400), [study PR #402](https://github.com/JeffreyEarly/wave-vortex-model/pull/402), and [advisory API PR #405](https://github.com/JeffreyEarly/wave-vortex-model/pull/405). InternalModes owns numerical assessment; these WVM-owned physical recipes remain outside the provider. This verification does not change production defaults.
