---
layout: default
title: IMDiscreteTransform
parent: IMDiscreteTransform
grand_parent: Discrete transforms
nav_order: 1
mathjax: true
---

#  IMDiscreteTransform

Create a scalar discrete transform from canonical matrices.


---

## Declaration
```matlab
 transform = IMDiscreteTransform(options)
```
## Parameters
+ `options.z`  physical sample points
+ `options.weights`  quadrature weights
+ `options.modeNumber`  retained mode labels
+ `options.normalization`  basis normalization name
+ `options.inverseMatrix`  inverse transform matrix containing the sampled modes
+ `options.metricMatrix`  sampled metric matrix
+ `options.targetGramMatrix`  continuous diagonal Gram target

## Returns
+ `transform`  initialized scalar discrete transform

## Discussion

  Let $$n_z$$ be the number of samples and $$n_m$$ the number of
  retained modes. `inverseMatrix` must be $$n_z\times n_m$$;
  `z` and `weights` must each contain $$n_z$$ entries;
  `metricMatrix` must be a symmetric $$n_z\times n_z$$ matrix;
  and `targetGramMatrix` must be a diagonal $$n_m\times n_m$$
  matrix with finite, nonzero diagonal entries. `modeNumber`
  supplies one label for each retained mode column.

  The constructor derives `gramMatrix`, `forwardMatrix`, and all
  quality diagnostics from these inputs. Most users should build
  this object from a solved basis with
  `IMBasisSet.discreteTransform`; direct construction is useful
  for alternative transform builders.

  ```matlab
  z = [-1; -0.5; 0];
  weights = [0.25; 0.5; 0.25];
  inverseMatrix = [1 0; 1 1; 0 1];
  metricMatrix = diag(weights);
  transform = IMDiscreteTransform(z=z,weights=weights,modeNumber=[1 2],normalization="unity",inverseMatrix=inverseMatrix,metricMatrix=metricMatrix,targetGramMatrix=eye(2));
  ```
