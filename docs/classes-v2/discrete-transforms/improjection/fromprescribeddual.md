---
layout: default
title: fromPrescribedDual
parent: IMProjection
grand_parent: Discrete transforms
nav_order: 7
mathjax: true
---

#  fromPrescribedDual

Construct a prescribed physical source dual without inventing synthesis data.


---

## Declaration
```matlab
 projection = IMProjection.fromPrescribedDual(samplePairingMatrix,sampleGram,targetGram,options)
```
## Parameters
+ `samplePairingMatrix`  finite real or complex column-by-sample pairing operator
+ `sampleGram`  real symmetric signed system for sampled coefficients
+ `targetGram`  real symmetric continuous coefficient system
+ `options.majorantGramMatrix`  positive coefficient error metric
+ `options.activeColumnMask`  active output coordinates, default all
+ `options.columnLabels`  distinct scientific output labels
+ `options.provenance`  required description of the supplied physical dual

## Returns
+ `projection`  source projection with no sampled synthesis basis

## Discussion

  `samplePairingMatrix` maps sampled source values into signed coefficient
  pairings. `sampleGram` supplies the signed coefficient system to solve;
  it need not equal a single observed component's Gram. `targetGram` is the
  corresponding continuous reference system. Prefixes restrict rows and
  coefficient systems together and solve each restricted system separately.

  The caller owns the physical recipe and supplies explicit provenance.
  The object can measure product errors using qualified reference pairings,
  but Gram and round-trip assessment are unavailable without a synthesis
  basis. Their diagnostics are NaN and capability flags are false.
