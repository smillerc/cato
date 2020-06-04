# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
