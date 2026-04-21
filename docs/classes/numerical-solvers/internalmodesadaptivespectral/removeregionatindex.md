---
layout: default
title: RemoveRegionAtIndex
parent: InternalModesAdaptiveSpectral
grand_parent: Classes
nav_order: 7
mathjax: true
---

#  RemoveRegionAtIndex

Merge neighboring adaptive regions by removing one boundary.

> Developer documentation: this item describes internal implementation details.


---

## Declaration
```matlab
 [newBoundaries,newSigns] = RemoveRegionAtIndex(oldBoundaries,oldSigns,index)
```
## Parameters
+ `oldBoundaries`  current region boundaries
+ `oldSigns`  signs of `N2-\omega^2` in each region
+ `index`  boundary index to remove

## Returns
+ `newBoundaries`  merged region boundaries
+ `newSigns`  updated regional signs

## Discussion

                  Given a list of boundaries (length n), and the signs of the
  regions created by the boundaries (length n-1), this removes
  the boundary at some index, and returns the appropriate signs
