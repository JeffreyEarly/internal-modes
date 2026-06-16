---
layout: default
title: modeFamily
parent: IMInternalModes
grand_parent: Core
nav_order: 14
mathjax: true
---

#  modeFamily

Physical mode-family declaration.


---

## Discussion

  `modeFamily` tells internal-mode utilities which physical
  catalog and coupled normalization rules are meaningful for this
  EVP. The default `"none"` installs only generic internal-mode
  behavior. The `"geostrophic"` family declares the hydrostatic
  geostrophic `F`/`G` family, enabling the generalized
  boundary-condition catalog and the coupled `geostrophic`
  normalization convention.
