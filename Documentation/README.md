Documentation
=============

+ The `WebsiteDocumentation` folder contains the hand-authored markdown source for the GitHub Pages site.
+ The `../tools/build_website_documentation.m` script rebuilds the published `../docs` folder by copying the website source, generating the version-history page from `../CHANGELOG.md`, rendering tutorials from `../Examples/Tutorials`, and extracting class reference pages with `class-docs`.
+ The generated `../docs` folder is committed output. Edit `Documentation/WebsiteDocumentation` and rerun the build script rather than editing `docs/` by hand.
