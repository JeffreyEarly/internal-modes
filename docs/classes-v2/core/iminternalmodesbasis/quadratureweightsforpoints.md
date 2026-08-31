---
layout: default
title: quadratureWeightsForPoints
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 16
mathjax: true
---

#  quadratureWeightsForPoints

Fit one quadrature rule to aligned internal-mode F/G channels.


---

## Declaration
```matlab
 [weights,weightFit] = quadratureWeightsForPoints(basisSet,options)
```
## Parameters
+ `options.z`  increasing fixed physical sample points
+ `options.nModes`  number of leading aligned family modes
+ `options.variables`  requested direct channels, F and/or G
+ `options.objective`  built-in stacked objective or custom callback
+ `options.nonnegative`  impose nonnegative shared weights
+ `options.constrainDepth`  impose full-depth shared weights

## Returns
+ `weights`  fitted shared quadrature weights
+ `weightFit`  specialized shared-fit diagnostics

## Discussion

Every requested directly representable channel contributes its compressed
normalized-Gram Frobenius system. The systems are stacked without an
additional channel weight and solved for one nonnegative, full-depth
constrained weight vector by default. Identically zero variable columns
are omitted from that variable's objective while their family positions
remain present in the returned transform.

A custom objective receives the stacked system in
`context.normalizedGramA` and `context.normalizedGramB`, row provenance in
`context.normalizedGramVariables`, and scalar-style per-variable contexts
in `context.variableContexts.F` and `.G`.
