# Analytical Validation Summary

This note summarizes the v2 analytical-solution validation tests in
`UnitTestsV2`. It is intentionally test-facing documentation, not generated API
documentation.

The analytical validation backbone currently lives in two suites:

| Test suite | Analytical basis | Main purpose |
| --- | --- | --- |
| `IMConstantStratificationValidationTests.m` | `IMBasisSet.constantStratification(...)` | Exact constant-stratification formulas, free-surface branches, solver comparisons, and normalization checks. |
| `IMExponentialStratificationValidationTests.m` | `IMBasisSet.exponentialStratification(...)` | Exponential Bessel-solution branches, free-surface branches, solver comparisons, and v1 oracle parity. |

## V1 Oracle Dependencies

These comparisons should vanish before v1 is removed. They currently validate
root values, equivalent depths, and mode shapes against the legacy analytical
classes. Replace them with self-contained formulas, root fixtures, or compact
precomputed reference data.

| Test suite | Test/helper | V1 class used | Purpose of comparison | Replacement needed when v1 is removed |
| --- | --- | --- | --- | --- |
| Constant | `freeSurfaceGModesMatchV1OracleAndBoundaryResidual` | `InternalModesConstantStratification` | Fixed-wavenumber free-surface `G` modes with `f0=0`; compares `h`, `G`, `F`, mode numbers, and boundary residuals. | Direct free-surface constant formulas for the surface branch and shifted interior roots. |
| Constant | `freeSurfaceWavenumberRegimesMatchV1Oracle` through `v1FreeSurfaceOracle` | `InternalModesConstantStratification` | Free-surface fixed-wavenumber regimes `0.1*kStar`, `kStar`, and `10*kStar`; compares `h`, surface-pressure-rescaled `G/F`, and boundary residuals. | Root fixtures or direct formulas for trig, transition, and hyperbolic free-surface branches. |
| Constant | `freeSurfaceFrequencyRegimesMatchV1SurfaceBranch` through `v1FreeSurfaceOracle` | `InternalModesConstantStratification` | Free-surface fixed-frequency regimes `0.1*N0`, `N0`, and `10*N0`; compares retained branches, surface-pressure-rescaled `G/F`, and boundary residuals. | Direct formulas or fixtures for low-frequency interior branches and high-frequency surface-only branches. |
| Exponential | `fixedWavenumberGModesMatchV1Oracle` through `v1ModesAtWavenumber` | `InternalModesExponentialStratification` | Rigid/rigid fixed-wavenumber `G` modes; compares `h`, `G`, `F`, mode numbers, and sign convention. | Self-contained Bessel-root fixtures for fixed-wavenumber rigid modes. |
| Exponential | `fixedFrequencyGModesMatchV1Oracle` through `v1ModesAtFrequency` | `InternalModesExponentialStratification` | Rigid/rigid fixed-frequency `G` modes below `N0`; compares `h`, `G`, `F`, mode numbers, and sign convention. | Self-contained Bessel-root fixtures for fixed-frequency rigid modes. |
| Exponential | `hydrostaticGModesMatchV1Oracle` through `v1ModesAtFrequency(..., omega=0)` | `InternalModesExponentialStratification` | Rigid/rigid hydrostatic `G` modes; compares `h`, `G`, `F`, sign convention, and rigid endpoint residuals. | Self-contained hydrostatic Bessel-root fixtures. |
| Exponential | `freeSurfaceWavenumberGModesMatchV1Oracle` through `v1FreeSurfaceModesAtWavenumber` | `InternalModesExponentialStratification` | Free-surface fixed-wavenumber `G` modes in three `kStar` regimes; compares `h`, surface-normalized `G/F`, mode numbers, and boundary residuals. | Free-surface Bessel-root fixtures for surface and interior branches. |
| Exponential | `freeSurfaceFrequencyGModesMatchV1Oracle` through `v1FreeSurfaceModesAtFrequency` | `InternalModesExponentialStratification` | Free-surface low-frequency and hydrostatic `G` modes; compares `h`, surface-normalized `G/F`, mode numbers, and boundary residuals. | Free-surface Bessel-root fixtures for fixed-frequency and hydrostatic limits. |
| Exponential | `highFrequencyFreeSurfaceGModesMatchV1Oracle` through `v1FreeSurfaceModesAtFrequency` | `InternalModesExponentialStratification` | Free-surface high-frequency `G` surface-only modes at `omega=N0` and `10*N0`; compares `h`, `G`, `F`, mode number `-1`, and boundary residuals. | Surface-branch Bessel-root fixtures for `omega >= N0`. |

## Constant-Stratification Coverage

| Scenario | EVP/factory | Boundaries | Normalization(s) | Compared quantities | Oracle/source | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| Rigid fixed-wavenumber `G` | `waveModesAtWavenumber(k=1e-4)` | surface `rigid`, bottom `rigid` | `wMax` | `h`, eigenvalues, `modeNumber`, `G(z)`, `F(z)` | Local direct formula helper | Checks mode numbers `1:nModes` and exact trigonometric shapes. |
| Rigid fixed-frequency `G` | `waveModesAtFrequency(omega=0.8*N0)` | surface `rigid`, bottom `rigid` | `wMax` | `h`, eigenvalues, `modeNumber`, `G(z)`, `F(z)` | Local direct formula helper | Covers fixed-frequency equivalent-depth relation below `N0`. |
| Rigid hydrostatic `G` | `hydrostaticGModes()` | surface `rigid`, bottom `rigid` | `wMax` | `h`, eigenvalues, `modeNumber`, `G(z)`, `F(z)` | Same direct formula as fixed-frequency with `omega=0` | Validates the hydrostatic `G` factory as the zero-frequency limit. |
| Rigid hydrostatic `F` | `hydrostaticFModes()` | surface `rigid`, bottom `rigid` | `uMax` | null mode, `h`, eigenvalues, `modeNumber`, `G(z)`, `F(z)`, diagnostic relation | Local direct `F` formula helper | Checks `modeNumber=0:(nModes-1)`, `h_0=Inf`, `lambda_0=0`, `F_0=1`, and `G_0=0`. |
| Rigid surface-pressure factors | all four rigid factories | surface `rigid`, bottom `rigid` | `surfacePressure` | normalization factors and `F(surface)=1` | Direct analytical factor formulas | Includes fixed-wavenumber `G`, fixed-frequency `G`, hydrostatic `G`, and hydrostatic `F`. |
| Hydrostatic `F` geostrophic Gram metrics | `hydrostaticFModes()` | surface `rigid`, bottom `rigid` | `geostrophic` | `gramMatrix("G")`, `gramMatrix("F")`, diagonals, off-diagonals | Exact constant-stratification metrics | Verifies null-mode metrics and baroclinic geostrophic normalization. |
| Free-surface fixed-wavenumber `G` | `waveModesAtWavenumber(..., surfaceBoundary=free)` | surface `free`, bottom `rigid` | default and `surfacePressure` depending on test | `h`, `G(z)`, `F(z)`, `modeNumber`, endpoint residuals | v1 oracle | Covers a representative `k=1e-4` case plus `0.1*kStar`, `kStar`, and `10*kStar`. |
| Free-surface fixed-frequency `G` | `waveModesAtFrequency(..., surfaceBoundary=free)` | surface `free`, bottom `rigid` | `surfacePressure` | retained `h`, `G(z)`, `F(z)`, `modeNumber`, endpoint residuals | v1 oracle | Covers `omega < N0`, `omega = N0`, and `omega > N0`; high-frequency cases retain only the surface mode. |
| Free-surface hydrostatic `G` | `hydrostaticGModes(surfaceBoundary=free)` | surface `free`, bottom `rigid` | `surfacePressure` | `h`, `G(z)`, `F(z)`, `modeNumber`, endpoint residuals | v2 fixed-frequency basis at `omega=0` | Confirms the hydrostatic free-surface branch is the fixed-frequency limit. |
| Rigid high-frequency rejection | `waveModesAtFrequency(omega>=N0)` | surface `rigid`, bottom `rigid` | not applicable | error id | Expected v2 unsupported error | Ensures degenerate rigid-surface interior modes are not silently returned. |
| Unsupported hydrostatic `F` boundaries | `hydrostaticFModes(surfaceBoundary=noSlip)` | surface `noSlip`, bottom `rigid` | not applicable | error id | Expected v2 unsupported error | Keeps constant analytical `F` support limited to rigid/rigid endpoints. |

## Exponential-Stratification Coverage

| Scenario | EVP/factory | Boundaries | Normalization(s) | Compared quantities | Oracle/source | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| Factory smoke | `hydrostaticGModes()` | surface `rigid`, bottom `rigid` | default | class, `N0`, `b`, `modeNumber`, `N2(surface)`, sign convention | v2 basis object | Confirms the factory returns `IMBasisSetExponentialStratification`. |
| Rigid fixed-wavenumber `G` | `waveModesAtWavenumber(k=1e-4)` | surface `rigid`, bottom `rigid` | `uMax` | `h`, `modeNumber`, `G(z)`, `F(z)`, `F(surface)>0` | v1 oracle | Tests ordinary Bessel roots and shapes. |
| Rigid fixed-frequency `G` | `waveModesAtFrequency(omega=0.8*N0)` | surface `rigid`, bottom `rigid` | `uMax` | `h`, `modeNumber`, `G(z)`, `F(z)`, `F(surface)>0` | v1 oracle | Tests fixed-frequency Bessel roots below `N0`. |
| Rigid hydrostatic `G` | `hydrostaticGModes()` | surface `rigid`, bottom `rigid` | `uMax` | `h`, `modeNumber`, `G(z)`, `F(z)`, sign convention, rigid residuals | v1 oracle at `omega=0` | Confirms surface and bottom `G` traces vanish. |
| Free-surface fixed-wavenumber `G` | `waveModesAtWavenumber(..., surfaceBoundary=free)` | surface `free`, bottom `rigid` | `uMax` with surface-normalized comparison | `h`, `modeNumber`, `G(z)`, `F(z)`, endpoint residuals | v1 oracle | Covers `0.1*kStar`, `kStar`, and `10*kStar`; expected labels are `[-1 1:(nModes-1)]`. |
| Free-surface low-frequency and hydrostatic `G` | `waveModesAtFrequency(..., surfaceBoundary=free)` and `hydrostaticGModes(surfaceBoundary=free)` | surface `free`, bottom `rigid` | `uMax` with surface-normalized comparison | `h`, `modeNumber`, `G(z)`, `F(z)`, endpoint residuals | v1 oracle | Covers `omega=0.1*N0` and `omega=0`. |
| Free-surface high-frequency `G` | `waveModesAtFrequency(..., surfaceBoundary=free)` | surface `free`, bottom `rigid` | `surfacePressure` | `h`, surface-only `modeNumber=-1`, `G(z)`, `F(z)`, endpoint residuals | v1 oracle | Covers `omega=N0` and `10*N0`. |
| Rigid hydrostatic `F` | `hydrostaticFModes()` | surface `rigid`, bottom `rigid` | `geostrophic`, `unity`, `surfacePressure` | null mode, baroclinic `h`, `G/F` shapes, endpoint residuals, `gramMatrix("F")`, `F(surface)=1` | v2 hydrostatic `G` baroclinic basis plus exact null-mode checks | Checks `modeNumber=0:(nModes-1)`, `lambda_0=0`, `h_0=Inf`, `F_0=1`, and `G_0=0`. |
| Single free-surface mode | `hydrostaticGModes(surfaceBoundary=free)` | surface `free`, bottom `rigid` | default | `modeNumber`, `h`, endpoint residuals | v2 analytical basis | Ensures `nModes=1` returns only the surface branch. |
| Surface-pressure normalization | fixed-wavenumber rigid and free-surface `G` | surface `rigid` or `free`, bottom `rigid` | `surfacePressure` | `F(surface)=1`; for free surface, `G(surface)=1`; for rigid surface, `G(surface)=0` | v2 analytical basis | Checks normalization behavior without v1. |
| Unsupported analytical boundaries | several factories | free-surface `F`, rigid too-fast frequency, no-slip bottom, free bottom, active surface | not applicable | error ids | Expected v2 unsupported errors | Keeps exponential analytical support limited to ordinary non-active `F/G` branches with rigid bottom. |
| Deferred dependency and old-token scan | analytical implementation files | not applicable | not applicable | source text | v2 source scan | Verifies no runtime `chebfun`, `rigidLid`, or `freeSlip` tokens in v2 analytical code. |

## Solver Validation

The same four numerical solvers are compared against the analytical v2 bases.
The tests check equivalent depths, mode numbers, mode shapes, sign convention,
boundary residuals, and, where listed, unity Gram-matrix behavior.

| Solver | Constant cases | Constant tolerances | Exponential cases | Exponential tolerances | Gram/unity checks |
| --- | --- | --- | --- | --- | --- |
| `IMSolverSpectral` | rigid fixed-wavenumber `G`; free-surface high-frequency surface branch; hydrostatic `F` | rigid and hydrostatic `h`: `1e-8`; shape: `1e-3`; free-surface branch uses at least `1e-3` | rigid fixed-wavenumber `G`; free-surface fixed-wavenumber `G` | rigid `h`: `1e-8`, shape: `1e-7`; free-surface uses at least `3e-3` | Constant rigid `G` unity `gramMatrix("G")`; exponential rigid/free-surface `G` unity Gram checks. |
| `IMSolverWKBSpectral` | same constant cases | same as spectral for constant tests | same exponential cases | rigid `h`: `1e-5`, shape: `1e-5`; free-surface uses at least `3e-3` | Same Gram checks as spectral, with WKB tolerances. |
| `IMSolverDensitySpectral` | same constant cases | same as spectral for constant tests | same exponential cases | rigid `h`: `1e-4`, shape: `2e-3`; free-surface uses at least `3e-3` | Same Gram checks as spectral, with density-coordinate tolerances. |
| `IMSolverFiniteDifference` | same constant cases | rigid and hydrostatic `h`: `1e-5`; shape: `1e-3`; free-surface branch uses at least `1e-3` | same exponential cases | rigid `h`: `1e-5`, shape: `1e-3`; free-surface uses at least `3e-3` | Constant rigid `G` unity Gram tolerance is looser; exponential Gram tolerance is `1e-4`. |

## Unsupported Cases

| Area | Case | Expected result |
| --- | --- | --- |
| Constant analytical basis | rigid-surface fixed-frequency `omega >= N0` | `IMBasisSetConstantStratification:UnsupportedFrequency` |
| Constant analytical basis | hydrostatic `F` with non-rigid supported surface variant in the test | `IMBasisSetConstantStratification:UnsupportedBoundary` |
| Exponential analytical basis | hydrostatic `F` with free surface | `IMBasisSetExponentialStratification:UnsupportedBoundary` |
| Exponential analytical basis | rigid-surface fixed-frequency `omega >= N0` | `IMBasisSetExponentialStratification:UnsupportedFrequency` |
| Exponential analytical basis | no-slip bottom or free bottom | `IMBasisSetExponentialStratification:UnsupportedBoundary` |
| Exponential analytical basis | active surface boundary appended to the EVP | `IMBasisSetExponentialStratification:UnsupportedBoundary` |
| Deferred dependencies and old names | `chebfun`, `rigidLid`, and `freeSlip` in v2 analytical source | forbidden by source scan |

## Notes For Future Cleanup

- v1 oracle comparisons are currently useful because they cross-check Bessel and
  free-surface root branches, but they should be replaced before v1 removal.
- Constant free-surface v1 shape comparisons are surface-pressure rescaled in
  the v2 tests where needed, because v1's free-surface `kConstant`
  normalization is not the intended long-term oracle.
- The v2 analytical source itself is expected to remain free of runtime Chebfun
  dependence; the exponential test suite only adds a local Chebfun path so the
  legacy v1 oracle can run when that dependency exists.
