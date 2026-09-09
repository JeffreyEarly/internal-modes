---
layout: default
title: assess
parent: IMBasisCollection
grand_parent: Discrete transforms
nav_order: 2
mathjax: true
---

#  assess

Measure a fixed requested page through the common result contract.


---

## Declaration
```matlab
 assessment = assess(collection,z,weights,options)
```
## Discussion

  Products and leakage are explicit prepared recipes as documented by
  IMProjection.assess. This operation neither changes scientific modes nor
  creates an independent reference solve. Apply tolerances afterward with
  checkBasisAssessment(assessment,...).
