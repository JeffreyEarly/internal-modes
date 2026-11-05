# Shared-grid requalification with DD neutral labels

The exact provisional default at 65 points, with smoothing length 39.0625 m, recovers **151 physical amplitudes plus two DD neutral-label amplitudes** under the joint-recovery study's unchanged requirements. This was measured on 2026-11-02 after aligning neutral labels with the manuscript. The earlier archive tested neighboring lengths but did not contain the exact default; its previously quoted 39 m count is now directly supported.

This is a retained-subspace recovery result for one specified physical problem and sample-based operator policy. It includes an independently refined modal reference, boundary resolution, separate energy/enstrophy errors, finite-amplitude kernel budgets and sample-noise limits. It does not establish recovery from arbitrary unresolved input, nonlinear resolution, or production WVM integration.

## Paired grid comparison

Both budgets reran all 18 archived grids and added the exact provisional default, for 38 cases. All 36 archived cases reproduce exactly the same nodes, weights, mandatory boundary-shape errors, physical-family reference qualification/refinement evidence and boundary pass/fail decisions. Every previously qualified physical count is unchanged. The 65-point comparison has 17 independently validated maxima, including the added exact default; the other two grids fail mandatory boundary resolution.

| Grid at 65 points | Historical physical amplitudes | DD physical amplitudes | DD neutral amplitudes |
| --- | ---: | ---: | ---: |
| Piecewise depth | 127 | 127 | 2 |
| Physical depth | 133 | 133 | 2 |
| Smoothed, 31.0659 m | 151 | 151 | 2 |
| Smoothed, 36.9438 m | 151 | 151 | 2 |
| Exact default, 39.0625 m | Not previously measured | 151 | 2 |
| Smoothed, 43.9338 m | 152 | 152 | 2 |
| Smoothed, 73.8876 m | 151 | 151 | 2 |
| Smoothed, 147.7752 m | 144 | 144 | 2 |

The best tested count remains 152 at 43.9338 m. The exact default is one physical amplitude below that result and remains on the broad high-recovery region. No unexpected DD regression or reason to change the provisional grid rule was found.

At 33 points all 19 grids fail the unchanged 1% mandatory boundary-shape gate, so no retained-count search is accepted. The exact default, 78.125 m, has 3.92325% error. The best tested smoothed case, 73.8876 m, has 2.71302%; piecewise depth has 1.31473%. The smoothing-length response is not monotone under the fixed local reconstruction, and the exact default is not the best 33-point candidate. None of these tested choices reaches the required accuracy. This is a finite experiment, not a proof that every possible 33-point grid or reconstruction must fail.

## Physical problem and unchanged protocol

The column is 1000 m deep, with a neutral surface cap on `[-200,0]` m, a quintic transition on `[-400,-200]` m and constant buoyancy frequency `N=0.005 s^-1` below. Specifically, `t=clip((-z-200)/200,0,1)` and `N2=2.5e-5*(6*t^5-15*t^4+10*t^3)`. Latitude is 45 degrees, wavelength 10 km, and generalized-energy coefficients are `g0=0`, `gd=-0.01`. The horizontal wavenumbers are `k=(2*pi/10000)*cos(pi/6)`, `l=(2*pi/10000)*sin(pi/6)`, with `g=9.81 m/s^2` and rotation rate `7.2921e-5 s^-1`. All families share one grid and its fixed positive physical-depth weights. Sample reconstruction uses degree-five local polynomials.

The new labels are the manuscript's DD sine profiles on the cap, with `m_j=j*pi/200`, zero displacement at both ends and the neutral-side derivative traces. The first two labels are required. Relative to the historical DN cap, their exact potential-enstrophy norms increase by factors 4 and 16/9 for the same displacement amplitudes. Their functions and norms have changed; comparing neutral-tail gains does not compare identical physical sources.

The physical search retains both surface-wave signs, every nonpositive APV mode and both zero-APV boundary responses. Additional structures are consecutive prefixes of internal waves, mixed-layer waves and positive APV modes. Physical amplitude count is `5+2*nInternal+2*nMixed+nAPV`; neutral counts are reported separately. Different maximum-count allocations need not contain the same mixed-layer prefix.

Construction inventories contain 56 internal, 26 mixed-layer, 98 positive APV and eight neutral structures. Construction reference resolutions are `[336 352 368 384 400 416 432]`, with interval Gauss rules `[1025 1281]`. Independent validation adds two structures in each family, uses reference resolutions `[368 384 400 416 432 448 464]` and quadrature `[1281 1537]`. The actual accepted construction and validation resolutions are 352 and 384, respectively.

No acceptance budget was changed:

| Gate | Requirement |
| --- | --- |
| Reference / quadrature accuracy | `1e-7` / `1e-9` |
| Mandatory component H1 shapes / endpoint responses | Each at most 1% |
| Supported energy / enstrophy amplitude errors | Each at most 1% |
| Energy / enstrophy kernel errors | Separate original finite-amplitude absolute budgets |
| Sample rank / conditioning | Relative cutoff `1e-12`; condition number at most `1e10` |
| Energy / enstrophy sample-noise gains | Each at most `1e3` |

Projectors are the study's direct physical pairings followed by joint overlap correction. Same-inventory round trips are therefore construction identities. Recovery qualification additionally requires independently refined modal samples. A maximum-count recovery claim also requires closed family sample bounds; adjacent one-family additions are checked. Omitted-mode responses remain diagnostics and do not select the smoothing length or alter acceptance. Ordinary manuscript projections without this joint correction do not inherit these counts.

## Exact default at 65 points

| Measurement | Result |
| --- | ---: |
| Smoothing length | 39.0625 m = 2.5 mean grid spacings |
| Neutral-cap nodes | 7 |
| Physical / neutral amplitudes | 151 / 2 |
| Internal / mixed-layer / positive APV structures | 46 / 6 / 42 |
| Mandatory boundary-shape error | 0.405877% |
| Independently checked energy error | `3.97118e-7` |
| Independently checked enstrophy error | `6.31908e-9` |
| Energy / enstrophy noise gains | 266.706 / 632.657 |
| Energy / enstrophy kernel-budget ratios | 0 / `7.51161e-5` |

The five mandatory physical amplitudes are additional to these structure counts. The sampled matrix has rank 153 and scaled condition number approximately `5.96e6`. A different count allocation requires its own check: the historical common choice `[16 8 16]` includes eight mixed-layer structures and is not a subset of this six-mixed-layer winner. The neighboring archived grids already rejected its overlap correction; the newly measured exact default also rejects it with DD labels.

The rejected common choice exposes a pre-existing projector limitation. Its samples have full rank 71, scaled condition number approximately `3.796e6`, and diagnostic minimum exact-recovery energy/enstrophy noise gains about 179.1/11.35, both below the 1000 budget. The weak-dual overlap fails the declared numerical rank test despite that sampled information. The diagnostic fitted inverse establishes retained-subspace identifiability; it does not establish the independent physical projection or omitted-mode accuracy. This warrants follow-up in the common transform construction and should not be interpreted as a demonstrated information limit of the grid. It does not invalidate the separately accepted `[46 6 42]` result.

The subsequent [transform investigation](TransformArchitectureInvestigation.md) identifies this fixed-allocation failure as numerical scaling of the weak-test equations. Invertible left preconditioning recovers the same weak operator at the unchanged rank tolerance and passes the independent retained-recovery gates; an independent QR audit confirms the result. This refines the diagnosis without changing the original-policy maximum counts recorded here. It does not remove the measured unresolved-input contamination or establish a newly optimized bandwidth.

At the maximum, the next two omitted APV structures cause energy/enstrophy amplitude contamination of approximately 34.78/6.48; the next two omitted internal structures cause energy contamination approximately 37.72. Across all available omitted labels, combined energy/enstrophy contamination is approximately 268.50/21.31, with a false-enstrophy kernel-budget ratio approximately 79685. These are normalized diagnostic gains, not independently passing recovery errors. Large omitted-mode contamination therefore remains a material limitation of using maximum retained counts. Changes in the false neutral enstrophy response are expected when both the label functions and their norms change; physical sources with the same labels can be compared separately from the changed neutral tails.

## Reduced physical bands

An additional exact-default check keeps both labels, all five mandatory physical amplitudes, points, weights, source scales and acceptance budgets fixed. Fresh construction and independently finer validation inventories qualify both reduced choices. The largest passing choice is reused from the main sweep rather than rebuilt.

| Physical counts `[internal mixed APV]` | Physical + neutral amplitudes | Identical maximum-tail energy / enstrophy response | Next-two-per-family energy / enstrophy response |
| --- | ---: | ---: | ---: |
| `[46 6 42]` | 151 + 2 | 268.502 / 21.3112 | 37.7165 / 6.48103 |
| `[41 5 37]` (floor of 90%) | 134 + 2 | 2.15777 / 0.677034 | 1.22061 / 0.145588 |
| `[36 4 33]` (floor of 80%) | 118 + 2 | 2.15242 / 0.621619 | 0.557902 / 0.0830194 |

The identical-tail comparison uses exactly the same 134 source amplitudes outside the maximum retained set, with identical DD functions and norms, for all three rows. The adjacent-tail definition changes as the retained band shrinks and is reported separately. Including every available omitted source for each reduced choice gives energy/enstrophy responses 2.18412/0.741937 at 134 physical amplitudes and 2.18229/0.624413 at 118. Thus a 10% per-family margin sharply improves physical-tail contamination but does not provide a small-error aliasing guarantee. The identical-tail false-enstrophy kernel-budget ratio remains approximately 1753–1759 for the reduced choices.

Reducing physical counts leaves the two retained label functions unchanged. Their omitted-neutral response is correspondingly unchanged: approximately 0.250299 for all available neutral tails and 0.0830188 for their next two omitted labels. A physical-band margin does not by itself control unresolved neutral labels.

## Reproduction and provenance

The research code is in the sibling `literature/wave-vortex-model-coordinate-grid` repository. Its new `matlab/runNeutralLabelGridRequalification.m` driver reuses the two local historical archive configurations, reruns every archived grid at 65 and 33 points, and adds the exact default at each point count. It preserves the historical archives and records dirty provider/study patches, all MATLAB source hashes, manuscript hash/revision and before/after runtime-source stability.

The driver uses the local ignored archives `data/joint-grid-comparison-20261031T090020` and `data/joint-grid-comparison-20261031T101629` for paired comparisons. Fresh calculations can instead be reproduced with `runJointRecoveryGridComparison`, without historical references or cached mode objects, using the physics, inventories and resolutions above. The archived smoothing ratios relative to `L0=147.77523446738553` m are `[1/8 2^(-2.625) 1/sqrt(32) 2^(-2.25) 1/4 2^(-1.75) 1/sqrt(8) 1/2 1/sqrt(2) 1 sqrt(2) 2 4 8 16]`; add `2.5*1000/((Nz-1)*L0)` to measure the exact default. Include the piecewise-depth, raw-bandwidth and physical-depth controls, with `commonCounts=[16 8 16]`.

Local result directory: `literature/wave-vortex-model-coordinate-grid/data/neutral-dd-grid-requalification-20261102`. It contains separate per-budget comparisons, independently passing winners, omitted-mode diagnostics, exact-default native-provider grid comparisons and source provenance. These generated study archives remain local and ignored by Git.

The reduced-band evidence is in its `default-band-margin/` subdirectory. `matlab/requalifyDDDefaultBandMargin.m` reproduces that check using a completed 65-point comparison, with configurable source and fresh output directories. The executed temporary helper and the persisted helper have identical scientific computation; only portable paths, output options and batch thread setup differ. Their source hashes and the setup-only comparison are recorded in the local verification ledger. A completed comparison is required to reuse the maximum-count controls without repeating that gate.

The study's affected neutral-projection, physics and neutral-diagnostic tests were migrated to fixed endpoint labels. All 23 focused methods passed, including the unchanged independent-reference class; no numerical tolerance changed. Code Analyzer was clean on the three changed test files and both new drivers. An independent final MATLAB audit checked every maximum's closed family bounds and unchanged gates, all 51 adjacent one-family additions rejected on both construction and refined-reference checks, exact archive reproduction, and both provider default grids. Provider and study default nodes and weights agree exactly at both budgets. Runtime source hashes and the manuscript hash remained unchanged during the main run. Earlier provider verification is recorded in [neutral-label verification](NeutralLabelModesVerification.md).
