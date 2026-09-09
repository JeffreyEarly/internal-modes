# Transform and product controls for unchanged nearby-kappa bases

`assessNeighborKappaProducts` uses the #19 inventory API to compare the same continuous fields at the actual requested κ through two sampled projections: the ordinary requested-κ basis and the unchanged neighboring-κ basis. It does not replace the requested source fields with neighboring or independently solved source modes.

Two numerical controls are declared explicitly: GG→G uses the basis-owned wave G pairing, including its required endpoint terms; FG→F uses a plain L2 pairing. The latter is a mathematical observation-space control because the public API does not supply a scalar physical wave F projection. These controls do not constitute WVM source channels or valid Fourier-interaction inventories.

Each retained band is assessed independently. Three outputs share one fixed sparse inventory, budget reservation, and factor cache: ordinary coefficients versus the continuous requested-basis reference; unchanged-neighbor coefficients versus that same reference; and the direct difference between the two complete-band coefficient operators versus zero. The difference operator is rebuilt for every actual band; it is never obtained by slicing the inverse of a larger band. It uses the prescribed-dual API and makes no Gram or round-trip claim. Product maxima measure ordinary aliasing, reuse-induced coefficient difference, and total error separately. Their maxima are not added as an equality. Per-product triangle inequalities and exact-repeat zero differences are checked internally.

All three outputs use the same positive requested-coordinate coefficient metric. Original basis labels and normalization names are retained in results. Matching uses eigenvalue rank with the stored orientation; no fitted cross-Gram map, rotation, sign change, or rescaling is applied to reused coefficients. A diagnostic cross-Gram map and the portion of neighboring synthesis outside the requested span are reported separately, preventing small coefficient differences from being presented as complete field accuracy.

The primary continuous projection evaluates the original requested basis on 513 points. A 1025-point quadrature comparison and the independently solved requested-κ basis on 1025 points provide separate reference discrepancies. Source products remain the same requested scientific functions on all grids. The original experiment's reference guard is required, and a maximum relative reference discrepancy of 1e-4 is an explicit additional experimental guard. Inconclusive references remain inconclusive. The difference-dual zero reference does not independently qualify the basis; qualification comes from the ordinary/neighbor comparisons and the retained solve evidence.

The declared bounded selection is cases 1, 12, 17, 24, 29, 36, 52, and 96, bands 8 and 16, and 65 sampled coordinates. It includes exact-repeat controls and κ offsets across all four profile families and both surface conventions. Each case/band/channel reserves at most 2000 individual output products, including zeros, before output preparation. The original continuous bases and detailed assessment MAT files remain outside the repository.

The helper requires the tested #19 provider API on the MATLAB path. No provider classes are copied into the study. Reproduction uses `assessNeighborKappaProducts(studyDirectory,outputDirectory)` after loading the provider and adding the study folder. These results authorize no production approximation and do not resolve independent pressure/momentum residuals, arbitrary superpositions, or complete model qualification.

## Bounded results

All 32 declared controls completed in 17.97 seconds. Their inventories reserved and evaluated 10,464 individual output products; numerical assessment accounted for 5.22 seconds and the additional synthesis diagnostics for 0.21 seconds. Each output examined 58 of 64 possible input pairs at band 8 and 160 of 256 at band 16. These are fixed sparse cumulative-anchor selections, not all-pair or physical-interaction guarantees.

Fifteen controls passed the declared experimental reference guard; seventeen remained inconclusive. The exact-repeat controls had zero reuse-induced coefficient difference, and their total errors equaled their ordinary baseline errors. Every evaluated finite product satisfied the common-metric triangle check.

| Case and control | Band | Baseline aliasing | Reuse coefficient difference | Total error | Reference guard |
| --- | ---: | ---: | ---: | ---: | --- |
| 1: constant, rigid, GGtoG | 8 | 5.51188e-05 | 0 | 5.51188e-05 | passed |
| 12: constant, rigid, GGtoG | 8 | 5.51188e-05 | 4.88749e-12 | 5.51188e-05 | passed |
| 12: constant, rigid, FGtoF | 8 | 0.00816059 | 0.145553 | 0.145819 | passed |
| 24: constant, free, GGtoG | 16 | 0.000511548 | 0.00148525 | 0.0014848 | passed |
| 36: exponential, rigid, GGtoG | 8 | 0.000217574 | 0.145769 | 0.145775 | passed |
| 36: exponential, rigid, GGtoG | 16 | 0.00303791 | 0.145925 | 0.145928 | passed |
| 52: pycnocline, rigid, GGtoG | 8 | 1.69356e-05 | 2.65843e-05 | 2.65382e-05 | inconclusive |
| 96: doubleSharp, free, GGtoG | 16 | 0.0778993 | 0.0344914 | 0.0806757 | inconclusive |

The constant-stratification rigid-lid example illustrates the importance of keeping basis components distinct: at a 10% κ offset, G reuse changed coefficients by only about 4.9e-12, while the unchanged F basis in the explicitly L2 control changed coefficients by about 0.146. This tests unchanged matrices; it does not rule out a separately derived exact rescaling or other structured reuse method. In the exponential rigid-lid 10% offset case, G reuse also produced about 0.146 coefficient error. The G synthesis component outside the requested span was about 0.056 at band 8 and 0.028 at band 16, so a coordinate map alone would hide meaningful field differences.

The selected pycnocline and double-sharp cases remain inconclusive and cannot support an acceptance claim. Some smooth-profile L2 F controls also remain inconclusive because their 513-to-1025-point quadrature discrepancy exceeds 1e-4. Source fields deliberately remain the original requested solve on every grid: the independent output-dual comparison does not itself qualify independent source-field error. The original solve diagnostics remain a separate required experimental guard.

Compact measurements are committed in [product-summary.json](product-summary.json), with [source provenance](product-source-provenance.json). Detailed evidence is retained at `/private/tmp/im21-products/selected/`, including `product-summary.json`, `source-provenance.json`, and one MAT file per control. The retained full experiment is `/private/tmp/im21-study/full/`. Large MAT files are not committed. These local paths are output locations, not external dependencies of the provider API.
