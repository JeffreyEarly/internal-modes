# Zero-APV form and rotation verification

This record describes the earlier integration. For the current normalization, exact-count and projection contracts, see [ModeCreationAPI](ModeCreationAPI.md), [BasisNormalization](BasisNormalization.md) and [CommonScalarProjection](CommonScalarProjection.md). Historical test counts and retired API names below are preserved as evidence of that earlier state.

Historical verification of the incoming implementation before the capable-basis integration. Code examples and class names below describe that historical API; current usage is documented in [the integration plan](APISimplificationPlan.md). These results do not replace verification of the integrated implementation.

This increment implements #49 from `InternalModesEVP` at `73e5e69` (merged PR #64). The common basis owns the physical forms and changes of continuous coefficient coordinates; the existing shared engine owns numerical recovery. No new EVP family, observation-grid design, source adapter or WVM dependency is introduced.

Verification used MATLAB R2026b Update 1 on MACA64, with `matlab -batch` outside the local sandbox and one computational thread. Applicable package, MATLAB, caching, documentation and focused agent guides were followed.

## Independent controls and review

The common API and shared pencil repair had separate authors. An independent test author supplied eight common-basis methods and six coordinate-edge methods. The implementation authors reviewed each other's mathematics, and the coordinator reviewed the production code, tests, records, collection integration and documentation. Audit findings were corrected before the final integration run.

The **24 new test methods** passed their focused checks: ten pencil methods, eight common-basis methods and six edge methods. The pencil author also ran the existing numerical and analytical lower-level zero-APV classes: **37 of 37 methods passed** in that combined gate. A complex test-fixture Cholesky failure was resolved by computing the same positive majorant before unitary observation phases; signed forms and their tolerances were preserved.

Physical references include independent continuous kinetic/potential-energy integration, the free-surface contribution, endpoint-response traces and the manuscript's signed boundary-depth relation. Numerical and analytical backends cover free-surface and rigid-lid configurations with zero, one and two active endpoints. Profile derivatives, selected/reordered/empty maps, normalization restoration, repeated/chunked collection pages and already rotated lower-level sources are checked. Runtime instrumentation confirms that coordinate changes reuse the computed endpoint responses without another solve.

Recovery controls use complex observations outside the represented span and verify synthesis, signed forms, physical duals, recovery operators and coefficient scales under a dense coordinate change. All three shared policies are covered. The extended-span control includes an omitted latent direction and another source component outside the full span, preserving the policy and observation/fitting factors rather than replacing recovery with retained least squares.

## Numerical corrections

The former pencil helper could return a repeated signed eigenspace with nonorthogonal columns. For `L=R=[2 1;1 -2]`, the right Gram was `[1 .5;.5 -1]` despite a zero pencil residual. It could also reject `L=R=[0 1;1 0]` because the eigensolver's individual vectors were right-metric null, although the complete repeated eigenspace admits a signed orthonormal frame.

The shared helper now normalizes repeated eigenspaces through their restricted signed metric, reorthogonalizes close distinct roots and requires the final pencil residual, signed right Gram and coordinate condition checks to pass. Independent scaling of the input forms removes dimensional rank/null floors. Binary-exponent restoration preserves physical eigenvalue units and rejects unrepresentable results. Tests cover zero, repeated and near-repeated roots, signed metrics, scales from `1e-250` to `1e250`, overflow/underflow limits, singular/defective/nonreal pencils and empty pages.

The audit also found positive roundoff in an exactly null physical-energy direction. In a fully neutral rigid-lid column of unit depth, with f0=1, g=2, k=.3, zero boundary coefficients and an eight-point physical-z spectral solve, the canonical energy is proportional to `[1 -1;-1 1]`. Its null direction has zero F and nonzero constant G. A rotated self-pairing near `3.16e-30` was previously accepted with a normalization divisor near `1.78e-15`. A physical-majorant self-ratio alone could not detect this cancellation. The new per-column canonical arithmetic bound rejects that unresolved norm at relative tolerance `1e-12`, without a dimensional absolute floor, changing the pairing or removing the column. An independently quadrature-tuned signed null is also rejected.

## Final integration and static gates

The final full V2 suite passed **690 of 690 methods**, with no failures or incomplete methods. All 21 changed MATLAB files matched their SHA256 snapshots throughout the final run. Code Analyzer reported zero findings across all 21 changed MATLAB files. After a final bounded inherited-name metadata correction, only its three affected files were checked again; all remained clean. Temporary `IMModeBasis` API extraction produced 42 Markdown files, including all 11 new public method/property entries. No generated website files were changed.

Source Markdown fences, 30 local links and supported structured API tokens passed checks. There is no configured `docs:check` task, so source checks and temporary API extraction supply the documentation gates. Whitespace, repository scope, unchanged manifest/dependencies and unchanged website/release artifacts were checked. No released package snapshot was edited, and the original dirty manuscript worktree's tracked diff remains unchanged.

With the declared dependencies and authoring repository on the MATLAB path, the portable integration command is:

```matlab
suite = matlab.unittest.TestSuite.fromFolder("UnitTestsV2",IncludingSubfolders=true);
runner = matlab.unittest.TestRunner.withTextOutput();
results = runner.run(suite);
assert(all([results.Passed]) && ~any([results.Incomplete]));
```

## Evidence limits

These checks establish the supported zero-APV coordinate and form contract, including explicit failure of unresolved norms and unsupported pencils. They do not certify arbitrary-function solve convergence, universal observability on an external grid, a larger retained band, new physical source projection semantics or rotations across inactive sectors. Direct developer wrappers require the same declared stratification as their background; the factory paths construct both from the same physical specification. See [forms and rotations](BasisRotations.md) for the public API and [shared recovery](SharedRecovery.md) for preparation and qualification requirements.
