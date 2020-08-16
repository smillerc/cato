# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.5.0] - 2020-08-12

### Added

- New/Improved TVD framework (still WIP)
- Improved pycato dataset loading with dask + xarray

---
## [1.4.0] - 2020-07-31

### Added

- Added the AUSMPW+, M-AUSMPW+, SLAU, SLAU2, SD-SLAU, SD-SLAU2

### Changes

- Boundary zones (ghost layers) defaults to 2
- Contour files default to ghost layers off and 64bit floats

---
## [1.3.0] - 2020-07-08


### Added

- Added the AUSM+ family of solvers, including AUSM+-up and AUSM+-up all-speed

### Fixed

- Fixed the memory leak in GFortran versions due to polymorphic function results

### Internals

 - Major rework of the solver implementations to allow the user to select different flux solvers, e.g. FVLEG or AUSM+-up

---
## [1.2.0] - 2020-06-04

### Added

- New/Improved flux and slope limiters
- Round-off error checking
- Improved integration tests, i.e. `sod_1d`, `sedov_2d`, etc. (more to come)
- Shared memory parallelization via OpenMP
- Residual convergence history output via `residual_hist.csv`

### Fixed

- Fixed gradient reconstruction and removed extraneous reconstruction methods
- Fixed a major bug in the timestep calculation via the CFL number
