---
layout: default
title: objectiveName
parent: IMQuadratureWeightFit
grand_parent: Discrete transforms
nav_order: 13
mathjax: true
---

#  objectiveName

Name of the least-squares objective used to fit the weights.


---

## Discussion

  The default value is `"normalizedGramFrobenius"`, which identifies
  the aggregate normalized Gram objective described in the class
  overview. A custom objective may provide its own name in the
  returned specification struct.
