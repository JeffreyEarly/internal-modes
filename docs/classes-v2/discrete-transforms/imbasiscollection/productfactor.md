---
layout: default
title: productFactor
parent: IMBasisCollection
grand_parent: Discrete transforms
nav_order: 9
mathjax: true
---

#  productFactor

Bind one collection page, variable, derivative, and coefficient function.


---

## Parameters
+ `collection`  continuous collection snapshot
+ `options.id`  distinct factor identity
+ `options.page`  requested collection page
+ `options.variable`  available continuous variable
+ `options.derivativeOrder`  supported vertical derivative order
+ `options.countRole`  retained or fixed; endpoints require fixed
+ `options.coefficient`  function of z multiplying each column

## Returns
+ `factor`  immutable metadata and deferred continuous evaluator

## Discussion

  Independent reference collections must be supplied and matched explicitly
  by the caller through a custom factor evaluator. This factory evaluates
  the same scientific basis on every requested quadrature grid.
