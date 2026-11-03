Documentation
=============

+ The `WebsiteDocumentation` folder contains the hand-authored markdown source for the GitHub Pages site.
+ The `../tools/build_website_documentation.m` script rebuilds the published `../docs` folder by copying the website source, generating the version-history page from `../CHANGELOG.md`, rendering tutorials from `../Examples/Tutorials`, and extracting class reference pages with `class-docs`.
+ The generated `../docs` folder is committed output. Edit `Documentation/WebsiteDocumentation` and rerun the build script rather than editing `docs/` by hand.

+ [Common scalar physical projection and recovery](CommonScalarProjection.md): supported family/channel pairings, literal endpoint observations, explicit recovery policies and qualification limits.
+ [Prescribed vertical sampling](PrescribedVerticalSampling.md): evaluate point values, true cell averages and cell integrals on caller-supplied layouts, then prepare generic numerical recovery.
