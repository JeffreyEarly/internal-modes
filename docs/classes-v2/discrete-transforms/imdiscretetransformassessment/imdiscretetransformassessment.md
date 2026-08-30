---
layout: default
title: IMDiscreteTransformAssessment
parent: IMDiscreteTransformAssessment
grand_parent: Discrete transforms
nav_order: 1
mathjax: true
---

#  IMDiscreteTransformAssessment

Create a scalar discrete-transform assessment.


---

## Declaration
```matlab
 assessment = IMDiscreteTransformAssessment(options)
```
## Parameters
+ `options.transform`  retained production transform
+ `options.candidateTransform`  full candidate transform
+ `options.weightFit`  quadrature-weight fit or empty
+ `options.requestedPointCount`  requested point count
+ `options.prefixDiagnostics`  per-prefix diagnostic table
+ `options.gramPolicy`  Gram policy result
+ `options.leakagePolicy`  leakage policy result
+ `options.quadraticAliasingPolicy`  quadratic policy result
+ `options.limitingPolicy`  limiting policy name
+ `options.retentionReason`  readable retention reason

## Returns
+ `assessment`  initialized assessment

## Discussion
