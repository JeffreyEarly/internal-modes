---
layout: default
title: allProductsPlan
parent: IMProductInventory
grand_parent: Discrete transforms
nav_order: 2
mathjax: true
---

#  allProductsPlan

Reserve all supplied products for bounded independent validation.


---

## Parameters
+ `options.interactionIds`  independently chosen supplied interactions
+ `options.prefixCounts`  increasing requested counts
+ `options.productBudget`  explicit reservation including zeros

## Returns
+ `plan`  bounded validation reservation

## Discussion

  This validation control measures the supplied inventory only. It adds no
  physical interactions, universal coverage, or superposition guarantee.
