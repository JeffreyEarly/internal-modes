---
layout: default
title: parameters
parent: IMEigenvalueProblem
grand_parent: Core
nav_order: 17
mathjax: true
---

#  parameters

Additional coefficient parameters.


---

## Discussion

  Fields are copied into the coefficient context so custom
  coefficient functions can read named values without new public
  properties. `parameters` is not the whole coefficient context; it
  is merged into the context provided by the solver and EVP.
  Standard internal-mode factories add fields `f0`, `g`, and
  `formulation`; wave-mode factories also add `k` or `omega`.

  ```matlab
  evp = IMEigenvalueProblem( ...
      p=@(z,ctx) ctx.alpha*ones(size(z)), ...
      parameters=struct("alpha",2));
  ```
