---
layout: default
title: IMProductInventory
parent: IMProductInventory
grand_parent: Discrete transforms
nav_order: 1
mathjax: true
---

#  IMProductInventory

Validate an inventory without invoking any evaluator.


---

## Parameters
+ `factors`  cell row of deferred factors
+ `products`  ordered valid product-family table
+ `outputs`  cell row of deferred output projections

## Discussion

  Factors require id, family, labels, ordinals, frequencySigns,
  countRole (retained/fixed), evaluate(z,columns,referenceId),
  and provenance. Outputs require id, family, labels, ordinals,
  countRole, prepare(grids), and provenance. Product table columns
  are interactionId, channel, factorA, factorB, output; numeric
  references index the corresponding cell rows. Optional complex
  coefficient multiplies the entire product on every grid.
