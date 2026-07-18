---
layout: default
title: hFromEigenvalue
parent: IMInternalModes
grand_parent: Core
nav_order: 9
mathjax: true
---

#  hFromEigenvalue

Equivalent-depth conversion function.


---

## Discussion

`hFromEigenvalue` maps retained eigenvalues to equivalent depths
for internal-mode basis sets. The handle has signature
`h = hFromEigenvalue(lambda)`, so
$$h_j=\texttt{hFromEigenvalue}(\lambda_j).$$
