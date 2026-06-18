---
layout: default
title: IMHydrostaticBoundaryCondition
has_children: false
has_toc: false
mathjax: true
parent: Core
grand_parent: Class documentation V2
nav_order: 4
---

#  IMHydrostaticBoundaryCondition

Convert a hydrostatic `F`/`G` endpoint law to canonical form.


---

## Declaration

<div class="language-matlab highlighter-rouge"><div class="highlight"><pre class="highlight"><code>classdef IMHydrostaticBoundaryCondition</code></pre></div></div>

## Overview

`IMHydrostaticBoundaryCondition` stores one physical hydrostatic
endpoint law:

$$
g\left[a+\frac{c}{gh}\right]F(z_\ell)
=
\left[b+\frac{d}{gh}+egh\right]G(z_\ell).
$$

The coefficients `a`, `b`, `c`, `d`, and `e` belong to this
physical endpoint law. They are not the raw constructor coefficients
of `IMBoundaryCondition`.

The known hydrostatic endpoint rules give endpoint additions to the
physical bilinear forms

$$
\langle F_i,F_j\rangle_F
=
\int_{z_b}^{z_s}F_i(z)F_j(z)\,dz
+
\sum_\ell \Delta_F^\ell,
$$

and

$$
\langle G_i,G_j\rangle_G
=
\frac{1}{g}\int_{z_b}^{z_s}N^2(z)G_i(z)G_j(z)\,dz
+
\sum_\ell \Delta_G^\ell.
$$

The table lists the endpoint additions $$\Delta_F^\ell$$ and
$$\Delta_G^\ell$$ for supported rows:

| Endpoint law | $$\Delta_F^\ell$$ | $$\Delta_G^\ell$$ |
| --- | --- | --- |
| $$gaF=bG$$ | $$0$$ | $$\eta_\ell\frac{b}{ga}G^i_\ell G^j_\ell$$ |
| $$\frac{c}{h}F=bG$$ | $$-\eta_\ell\frac{c}{b}F^i_\ell F^j_\ell$$ | $$0$$ |
| $$g\left(a+\frac{c}{gh}\right)F=0$$ | $$0$$ | $$0$$ |
| $$g\left(a+\frac{c}{gh}\right)F=bG$$ | $$-\eta_\ell\frac{c}{b}F^i_\ell F^j_\ell$$ | $$\eta_\ell\frac{ga}{b}F^i_\ell F^j_\ell$$ |
| $$\frac{c}{h}F=\frac{d}{gh}G$$ | $$0$$ | $$\eta_\ell\frac{d}{gc}G^i_\ell G^j_\ell$$ |
| $$\left(b+\frac{d}{gh}\right)G=0$$ | $$0$$ | $$0$$ |
| $$(b+egh)G=0$$ | $$0$$ | $$0$$ |
| $$gaF=eghG$$ | $$-\eta_\ell\frac{a}{e}F^i_\ell F^j_\ell$$ | $$0$$ |
| $$\frac{c}{h}F=\left(b+\frac{d}{gh}\right)G$$ | $$-\eta_\ell\frac{1}{bc}\left(cF^i_\ell-\frac{d}{g}G^i_\ell\right)\left(cF^j_\ell-\frac{d}{g}G^j_\ell\right)$$ | $$\eta_\ell\frac{d}{gc}G^i_\ell G^j_\ell$$ |
| $$gaF=(b+egh)G$$ | $$-\eta_\ell\frac{1}{ae}\left(aF^i_\ell-\frac{b}{g}G^i_\ell\right)\left(aF^j_\ell-\frac{b}{g}G^j_\ell\right)$$ | $$\eta_\ell\frac{b}{ga}G^i_\ell G^j_\ell$$ |

Conversion to canonical EVP boundaries is handled by
`canonicalBoundary`. For an `F`-formulated EVP, the law must have
$$e=0$$ and converts as

```matlab
boundary = IMBoundaryCondition(a=-law.a,b=law.b,c=-law.c/g,d=law.d/g);
```

For a `G`-formulated EVP, the law must have $$d=0$$ and converts as

```matlab
boundary = IMBoundaryCondition(a=law.e,b=law.a,c=law.b/g,d=law.c/g);
```

If the requested formulation would put both $$h$$ and $$1/h$$ on the
same side of the endpoint law, the conversion is nonlinear in
$$\lambda=1/h$$ and cannot be represented by one canonical linear
boundary row.

After conversion, `IMInternalModes.innerProduct("F")` and
`IMInternalModes.innerProduct("G")` report which bilinear forms are
available for the resolved canonical boundary condition.

```matlab
g = 9.81;
law = IMHydrostaticBoundaryCondition(a=A/g,b=1);
surfaceBoundary = law.canonicalBoundary(formulation="G",g=g);
evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-4000 0],g=g,surfaceBoundary=surfaceBoundary);
```




## Topics
+ Create hydrostatic boundary laws
  + [`IMHydrostaticBoundaryCondition`](/internal-modes/classes-v2/core/imhydrostaticboundarycondition/imhydrostaticboundarycondition.html) Create a physical hydrostatic endpoint law.
+ Convert hydrostatic boundary laws
  + [`canonicalBoundary`](/internal-modes/classes-v2/core/imhydrostaticboundarycondition/canonicalboundary.html) Convert to a canonical scalar boundary condition.
+ Inspect hydrostatic boundary laws
  + [`a`](/internal-modes/classes-v2/core/imhydrostaticboundarycondition/a.html) Coefficient multiplying `F` on the \(h^0\) side.
  + [`b`](/internal-modes/classes-v2/core/imhydrostaticboundarycondition/b.html) Coefficient multiplying `G` on the \(h^0\) side.
  + [`c`](/internal-modes/classes-v2/core/imhydrostaticboundarycondition/c.html) Coefficient multiplying `F` on the \(1/h\) side.
  + [`d`](/internal-modes/classes-v2/core/imhydrostaticboundarycondition/d.html) Coefficient multiplying `G` on the \(1/h\) side.
  + [`e`](/internal-modes/classes-v2/core/imhydrostaticboundarycondition/e.html) Coefficient multiplying `G` on the \(h\) side.


---