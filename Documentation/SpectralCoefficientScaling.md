# Spectral coefficient scaling

Long external waves exposed a numerical limitation in the row-equilibrated Chebyshev pencil: eigenpairs with small normwise matrix residuals could still have inaccurate physical derivatives. At 100 km wavelength, 1000 m depth and exponential stratification, the existing WVM vertical-momentum check measured approximately `2.41e-7` at 64 coefficients against its `2e-7` allowance. Increasing the eigenproblem to 128 coefficients increased the residual to approximately `3.15e-6` on a 257-point evaluation grid. Direct differentiation of the native expansion reproduced the sampled-pressure residual, so the defect was not caused by the diagnostic's evaluation grid alone.

## Numerical change

For a canonical pencil `A*c = lambda*B*c`, retain the existing equation row equilibration `R` and solve instead

$$
(R A S)y = \lambda (R B S)y,\qquad c = S y,\qquad S_{jj} = \max(1,j)^{-2},
$$

where `j` is the zero-based Chebyshev polynomial degree. This invertible diagonal change of unknowns reduces high-degree coefficient roundoff that differentiation can amplify. It does not remove polynomial columns, alter the boundary equations, change retained counts or replace the resolved modes. The physical matrices remain available for classification and diagnostics.

`IMSolverSpectral` supplies the degree scaling through a protected solver hook. Other discretizations inherit unit scaling, preserving the finite-difference unknowns. The shared scalar/bulk solve recovers `c` before finite-generalized-mode filtering, sign orientation and normalization. Filtering uses the recovered coefficients with the original row-scaled metric matrix; it must not combine scaled eigenvectors with an unscaled metric, or vice versa.

## Verification and interpretation

`IMSpectralWaveMomentumTests` exercises all six requested wave modes at 100 km and 1 km wavelength, with constant and exponential stratification, `z`, WKB and density coordinates, and 64/128 eigenproblem coefficients. It independently differentiates native `G` twice to obtain `g*h*G_zz`, and separately differentiates sampled `F` using a 257-point physical derivative rule. It checks endpoint conditions and compares constant-stratification modes with the analytical solution. The direct residual allowance is `2e-8`; the sampled-field check preserves WVM's existing `2e-7` allowance. Eleven of the twelve parameterized cases fail on the unmodified beta.2 solver; all twelve pass after scaling.

The physical residual is the depth-weighted norm of `g*F_z + (N2-omega^2)*G`, divided by the sum of the separate pressure-gradient, inertia and buoyancy norms. This is a local equation-balance diagnostic, not a relative pressure-field error bound. Inferring `F_z` from that same equation would make the residual circular and is deliberately avoided.

The exploratory exponential-profile external-mode comparison gives direct residuals `1.70e-11` at 64 coefficients and `3.81e-11` at 128 after degree scaling. An independent physical-coordinate shooting solution agrees with the 64-point equivalent depth and pressure-mode shape to approximately `1.8e-13` and `2.6e-13`, respectively. These are bounded numerical comparisons, not universal guarantees or exact-solution error bounds.

The change does not introduce a public tuning parameter or automatic mode selection. Broader per-mode solve-quality diagnostics remain [InternalModes #10](https://github.com/JeffreyEarly/internal-modes/issues/10); a small matrix residual alone remains insufficient evidence of physical accuracy. The historical [API adoption verification](APIAdoptionVerification.md) records the beta.2 baseline limitation before this correction.

## Integration checks

On MATLAB R2025b Update 4, all **360 InternalModes V2 tests pass**, including the 12 new parameterized tests and existing canonical/analytical, APV, mean-density, negative/zero-mode, finite-difference, projection, and scalar/bulk regressions. All **76 selected WVM tests pass** on merged v5 revision `f67573ac`, covering the original free-surface wave qualification, exact bulk construction, mixed Boussinesq transform/evolution, resolution transfer/restart, and QG transform/calculus/diagnostics/conservation. The previously failing `variableModesAgreeWithIndependentShooting` test passes without changing its tolerance.

The downstream authoring advisory suite passes 17 of 18 cases. Its only failing assertion compares newly solved modes with numerical error estimates frozen under the old provider. The first three constant-profile prefix errors improve from `[7.9853e-7, 8.7664e-7, 1.1294e-6]` to `[1.0257e-7, 1.0257e-7, 1.9841e-7]`, exceeding that historical equality check's `1e-7` allowance. The five higher-prefix errors retain the old values within that allowance; count recommendations, explicit rejection, independent Gram/quadratic decisions, reference gating and budget checks all pass. Preserve the historical fixture and refresh the applicable numerical expectation with explicit provider provenance during the next WVM dependency adoption. This change does not relax the equality tolerance or rewrite the earlier study's results.

These WVM checks use the authoring provider correction; they are not evidence of a released or clean-installed new dependency. No package version, dependency floor, or released snapshot changes in this PR. The provider must be released/exported before WVM adopts it. Documentation was generated once, with only the intended version-history change. Code Analyzer introduces no new findings; the spectral solver retains its pre-existing unused `zBounds` argument diagnostic.
