---
layout: default
title: IMProjectionRecipe
parent: IMProjectionRecipe
grand_parent: Discrete transforms
nav_order: 1
mathjax: true
---

#  IMProjectionRecipe

Store a basis-owned recipe without changing its scientific state.


---

## Parameters
+ `options.available`  continuous pairing availability
+ `options.reason`  unavailable pairing explanation
+ `options.variable`  sampled scalar variable
+ `options.columnLabels`  scientific labels in column order
+ `options.normalization`  frozen normalization convention
+ `options.targetGramMatrix`  continuous signed Gram target
+ `options.majorantGramMatrix`  positive coefficient metric
+ `options.supportsLeakage`  family leakage capability
+ `options.supportsQuadratic`  family quadratic capability
+ `options.provenance`  representation and reference quality
+ `options.zDomain`  physical bounds
+ `options.evaluateFunction`  normalized basis evaluator
+ `options.weightFunction`  bound interior-weight evaluator
+ `options.spec`  basis-owned endpoint pairing specification

## Returns
+ `self`  immutable recipe
