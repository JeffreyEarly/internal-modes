---
layout: default
title: Inner products and norms
parent: Class documentation V2
nav_order: 4
permalink: /classes-v2/inner-products-and-norms
mathjax: true
---

# Inner products and norms

Internal-mode families with eigenvalue-dependent endpoint conditions live naturally in a Pontryagin space. Their scientific pairing can therefore be signed even though the interior weight is positive. InternalModes keeps that signed pairing for modal algebra and exposes its induced positive Hilbert majorant for numerical magnitudes.

For endpoint functionals $$L_\ell$$ and real endpoint coefficients $$\alpha_\ell$$, the signed pairing is

$$
[u,v]=\int_{z_b}^{z_s}w(z)\,\overline{u(z)}v(z)\,dz+\sum_\ell \alpha_\ell\overline{L_\ell[u]}L_\ell[v].
$$

The induced Hilbert majorant keeps the positive interior term and changes only the endpoint coefficients:

$$
(u,v)_+=\int_{z_b}^{z_s}w(z)\,\overline{u(z)}v(z)\,dz+\sum_\ell |\alpha_\ell|\overline{L_\ell[u]}L_\ell[v].
$$

This construction applies to both solved-variable boundary functionals and catalog endpoint terms. A family can have up to two negative endpoint directions. If every endpoint coefficient is nonnegative, the signed and majorant forms coincide.

## Which form to use

| Task | Use | APIs |
| --- | --- | --- |
| Modal orthogonality and signatures | Signed pairing | `innerProduct`, `gramMatrix`, `spectrum` |
| Projection functionals and coefficient recovery | Signed pairing | `modeProjectionFunctional`, `transformForward`, `targetGramMatrix` |
| Signed physical invariants | Signed pairing | family-specific Gram matrices and `spectrum` |
| Error tolerances and convergence tests | Hilbert majorant | `majorantGramMatrix`, `majorantNorm` |
| Aliasing magnitudes and state-size comparisons | Hilbert majorant | `targetMajorantGramMatrix`, quadratic assessment diagnostics |

Do not use $$\sqrt{[u,u]}$$ or $$\sqrt{|[u,u]|}$$ as a norm for an arbitrary state. Positive and negative directions can cancel in $$[u,u]$$, and the absolute value does not restore the triangle inequality. The property `innerProductNormFactor` is a per-mode normalization convention, $$\sqrt{|[V_j,V_j]|}$$; it is not an arbitrary-state norm.

## Basis APIs

Numerical `IMInternalModesBasis` and analytical `IMAnalyticalInternalModesBasis` objects provide parallel interfaces:

```matlab
signedRecipe = basis.evp.innerProduct("G");
majorantRecipe = basis.majorantInnerProduct(variable="G");

Gamma = basis.gramMatrix(variable="G");
Mplus = basis.majorantGramMatrix(variable="G");
stateSize = basis.majorantNorm(coefficients,variable="G");
```

Both Gram methods accept `zBounds` when a partial-depth matrix is required. `majorantNorm` evaluates

$$
\lVert u\rVert_+=\sqrt{\mathbf c^*M_+\mathbf c}.
$$

The modal basis diagonalizes the signed invariant in the usual normalization, but it need not diagonalize the majorant. Consequently there is no “majorant spectrum.” `spectrum` remains the signed modal invariant and may contain negative entries.

## Discrete transforms and quadratic certification

`IMInternalModesDiscreteTransform.targetGramMatrix(variable=...)` stores the continuous signed target used by the projection contract. `targetMajorantGramMatrix(variable=...)` stores the corresponding positive target used to measure errors.

For an APV product $$p_{ij}$$, coupled quadratic certification evaluates

$$
q_N=\max_{i,j}\frac{\left\lVert\Pi_N^\mathrm{discrete}p_{ij}-\Pi_N^\mathrm{continuous}p_{ij}\right\rVert_+}{\lVert p_{ij}\rVert_+}.
$$

Both projections use signed pairings and signed Gram solves. Only the numerator and denominator magnitudes use the majorant. Negative APV modes keep their negative `modeNumber` and participate in every applicable `FF->F`, `FG->G`, and `GG->F` product pair. Assessments record `quadraticAliasingPolicy.projectionPairing="signedPontryagin"` and `quadraticAliasingPolicy.errorNorm="inducedHilbertMajorant"`.

Signed rejected-mode leakage remains unsupported because it does not yet have a separate scientific error contract. Coupled MDA quadratic products also remain unsupported.
