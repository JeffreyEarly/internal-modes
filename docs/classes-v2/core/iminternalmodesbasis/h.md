---
layout: default
title: h
parent: IMInternalModesBasis
grand_parent: Core
nav_order: 11
mathjax: true
---

#  h

Equivalent depths for the retained internal modes.


---

## Discussion

These are computed from the parent EVP as
$$h_j=\texttt{evp.hFromEigenvalue}(\lambda_j),$$
where $$\lambda_j$$ is `eigenvalues(j)`.
