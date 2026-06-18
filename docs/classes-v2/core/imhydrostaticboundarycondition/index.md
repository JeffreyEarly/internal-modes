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
endpoint law,
$$g\left[a+\frac{c}{gh}\right]F(z_\ell)
=
\left[b+\frac{d}{gh}+egh\right]G(z_\ell).$$
It is a helper for creating the canonical `IMBoundaryCondition`
required by `IMInternalModes` factories. EVPs still store only
`IMBoundaryCondition` objects.

For an `F`-formulated EVP, the law must have $$e=0$$ and converts to
$$\left(a_\mathrm{code},b_\mathrm{code},c_\mathrm{code},d_\mathrm{code}\right)
=
\left(-a,b,-\frac{c}{g},\frac{d}{g}\right).$$
For a `G`-formulated EVP, the law must have $$d=0$$ and converts to
$$\left(a_\mathrm{code},b_\mathrm{code},c_\mathrm{code},d_\mathrm{code}\right)
=
\left(e,a,\frac{b}{g},\frac{c}{g}\right).$$
If the requested formulation would put both $$h$$ and $$1/h$$ on the
same side of the endpoint law, the conversion is nonlinear in
$$\lambda=1/h$$ and cannot be represented by one canonical linear
boundary row.

Known physical endpoint additions from the hydrostatic catalog:

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