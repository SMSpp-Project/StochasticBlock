# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

### Changed

### Fixed

## [0.6.0] - 2025-12-12

### Added

- makefile-s and -c

- comprehensive test suite for DiscreteScenarioSet scenario reduction

- ScenarioReductionConfig class following SMS++ Configuration pattern

- [big] implemented scenario reduction functionality in
  DiscreteScenarioSet using CapacitatedFacilityLocationBlock

- unit tests of DiscreteScenarioSet

- [huge] DiscreteScenarioSet implementation of ScenarioGenerator

- [huge] new general ScenarioGenerator interface

### Changed

- optimized compute\_transport\_cost\_matrix (only computes upper
  triangle coefficients)

- better way to do sampling without replacement using
  std::discrete\_distribution which does sampling with replacement
  and keeping new indices in an unordered set

- wrapped debug couts with NDEBUG flag, refactoring tests

- StochasticBlock now depends on CapacitatedFacilityLocationBlock to
  implement scenario reduction capabilities

- adapted to new standard organization of makefiles

### Fixed

- BlockConfig ownership issue in set\_scenario\_reduction\_config

- compute\_scenario\_distance to properly compute ell Wasserstein
  of the euclidean norm by default.

- makefiles

## [0.5.0] - 2024-02-29

### Added

### Changed

### Fixed

## [0.4.3] - 2024-02-29

### Changed

- Adapt to new CMake / makefile organisation

## [0.4.2] - 2022-07-01

### Added

- Define the sense of the Objective of the StochasticBlock.

### Changed

- Default argument for "issueAMod" parameter in set_data() becomes eNoBlck

- Update load() and print() interfaces.

## [0.4.1] - 2021-12-08

### Added

- Makefile.

### Changed

- CMake file.

## [0.4.0] - 2021-05-02

### Changed

- Maintenance release.

## [0.3.0] - 2020-09-16

### Changed

- Maintenance release.

## [0.2.0] - 2020-03-06

### Added

- (De)serialization considering a vector of SimpleDataMappingBase.

## [0.1.0] - 2020-01-02

### Added

- First test release.

[Unreleased]: https://gitlab.com/smspp/stochasticblock/-/compare/0.6.0...develop
[0.6.0]: https://gitlab.com/smspp/stochasticblock/-/compare/0.5.0...0.6.0
[0.5.0]: https://gitlab.com/smspp/stochasticblock/-/compare/0.4.3...0.5.0
[0.4.3]: https://gitlab.com/smspp/stochasticblock/-/compare/0.4.2...0.4.3
[0.4.2]: https://gitlab.com/smspp/stochasticblock/-/compare/0.4.1...0.4.2
[0.4.1]: https://gitlab.com/smspp/stochasticblock/-/compare/0.4.0...0.4.1
[0.4.0]: https://gitlab.com/smspp/stochasticblock/-/compare/0.3.0...0.4.0
[0.3.0]: https://gitlab.com/smspp/stochasticblock/-/compare/0.2.0...0.3.0
[0.2.0]: https://gitlab.com/smspp/stochasticblock/-/compare/0.1.0...0.2.0
[0.1.0]: https://gitlab.com/smspp/stochasticblock/-/tags/0.1.0
