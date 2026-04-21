---
layout: default
title: Installation
nav_order: 2
description: Installation instructions
permalink: /installation
---

# Installation

Install `Internal Modes` either through the OceanKit MPM repository or directly from the authoring repository.

## Runtime requirements

- MATLAB R2024b or newer
- the `SplineCore`, `Distributions`, and `chebfun` dependencies used by the package implementations

## Install from OceanKit

Clone the OceanKit repository:

```text
git clone https://github.com/JeffreyEarly/OceanKit.git
```

Then register it with MPM from within MATLAB:

```matlab
mpmAddRepository("OceanKit", "path/to/OceanKit")
mpminstall("InternalModes")
```

## Install from the authoring repository

Clone the source repository:

```text
git clone https://github.com/JeffreyEarly/internal-modes.git
```

Then install it from MATLAB:

```matlab
mpminstall("local/path/to/internal-modes", Authoring=true)
```

## Development and documentation

If you want to rebuild the website documentation locally, make sure the OceanKit `ClassDocumentation` or `class-docs` tooling is on the MATLAB path in addition to the runtime dependencies.
