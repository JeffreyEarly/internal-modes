---
layout: default
title: metadata
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 14
mathjax: true
---

#  metadata

Additional EVP metadata.


---

## Discussion

  Fields are copied into the coefficient context so custom
  coefficient functions can read scalar parameters without new
  public properties. Standard internal-mode factories add fields
  `f0`, `g`, and `formulation`; wave-mode factories also add `k` or
  `omega`. User-supplied metadata fields are preserved on
  `evp.metadata` and are visible to coefficient handles through
  `ctx.fieldName`.
