---
layout: default
title: withDiagnostics
parent: InternalModesTransform
grand_parent: Classes
nav_order: 47
mathjax: true
---

#  withDiagnostics

Return a copy with updated selection diagnostics.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 transform = withDiagnostics(self,options)
```
## Parameters
+ `self`  InternalModesTransform instance
+ `options.selectionReason`  reason retained modes were selected
+ `options.projectionResolvedModes`  projection-quality mode count
+ `options.nonlinearAliasLimit`  nonlinear aliasing mode limit

## Returns
+ `transform`  copied transform with diagnostics updated

## Discussion
