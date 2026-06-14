---
layout: default
title: h
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 4
mathjax: true
---

#  h

Equivalent depths for the retained internal modes.


---

## Discussion

  For numerical internal-mode solves, these are computed from the
  parent EVP as
  $$h_j=\texttt{evp.hFromEigenvalue}(\lambda_j),$$
  where $$\lambda_j$$ is `eigenvalues(j)`. Analytical basis classes
  may pass exact equivalent depths directly.
