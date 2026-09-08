---
layout: default
title: evaluate
parent: IMBasisCollection
grand_parent: Discrete transforms
nav_order: 5
mathjax: true
---

#  evaluate

continuous structures on selected requested pages.


---

## Parameters
+ `z`  physical coordinate column
+ `options.variable`  scalar variable name
+ `options.derivativeOrder`  requested vertical derivative order
+ `options.pages`  requested-page indices, empty means all
+ `options.columns`  explicit scientific column indices

## Returns
+ `values`  evaluated continuous basis columns by requested page

## Discussion

  The output is `nZ x nColumns x nSelectedPages` (MATLAB omits trailing
  singleton dimensions). Empty `pages` and `columns` select all entries.
  An omitted variable selects the first variable of the first selected basis.
  Explicit page and column order, including repeated indices, is preserved.
  Endpoint columns may be selected explicitly; they are never treated as an
  automatically accepted or truncated modal prefix.

  First derivatives use the existing solved-variable derivative evaluator.
  Boundary F derivatives use $$F_z=-N^2G/g$$. Other diagnostic derivatives
  require a separately implemented mathematical relation and are rejected.
