# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

### Changed

- the test of the discrete distribution removes the files it leaves behind
  with `std::filesystem`, rather than with a call to the shell that says
  nothing when it fails
- whoever links the module keeps it: the classes of a module register
  themselves in the factory from a static initialiser, and a linker that
  drops what looks unused takes the registration away with it, so the target
  now tells whoever links it to keep the symbol that forces the module in,
  and on ELF, where naming the symbol is not enough, the library as a whole
### Fixed

- the module does not depend on CapacitatedFacilityLocationBlock: the
  sub-problem of the scenario reduction is built through the factory, so
  neither the makefiles nor the packages require that module

## [0.7.0] - 2026-09-12

### Added

- `ScenarioReductionBlock`, and `DiscreteScenarioSet::load_from_memory()` and
  `get_scenario()`: a DiscreteScenarioSet is a container of the data, whose
  reduction is left to a Solver such as those of ScenarioReductionSolver

- `IndependentMultiStageScenarioGenerator`, the multi-stage scenario generator
  whose stages are independent

### Changed

- StochasticBlock no longer links CapacitatedFacilityLocationBlock, which
  nothing in it uses since the scenario reduction moved to
  ScenarioReductionSolver

- a multi-stage scenario generator is navigated through Views only: the
  internal cursor is gone, `descend()` and `climb()` move the View itself,
  which has `clone()`

- (re-)defining the pool is the only write on a scenario generator: a View
  requires an initialized pool and detects having been invalidated by a new
  one

- the pool of a node is reduced through the DiscreteScenarioSet machinery, the
  node keeping its children as the universe and only restricting the pool the
  Views read; this is a local heuristic, not the selection of scenarios a tree
  asks for

- the version of the module is the git tag of its repository, or the
  VERSION.txt of a release tarball, and the shared library carries it: its
  SONAME is major.minor while the major is 0, and it is installed with an
  RPATH relative to itself, so that an installed tree keeps working wherever
  it is moved

### Fixed

- `ScenarioGenerator::serialize()` writes the type attribute

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

[Unreleased]: https://gitlab.com/smspp/stochasticblock/-/compare/0.7.0...develop
[0.7.0]: https://gitlab.com/smspp/stochasticblock/-/compare/0.6.0...0.7.0
[0.6.0]: https://gitlab.com/smspp/stochasticblock/-/compare/0.4.3...0.6.0
[0.4.3]: https://gitlab.com/smspp/stochasticblock/-/compare/0.4.2...0.4.3
[0.4.2]: https://gitlab.com/smspp/stochasticblock/-/compare/0.4.1...0.4.2
[0.4.1]: https://gitlab.com/smspp/stochasticblock/-/compare/0.4.0...0.4.1
[0.4.0]: https://gitlab.com/smspp/stochasticblock/-/compare/0.3.0...0.4.0
[0.3.0]: https://gitlab.com/smspp/stochasticblock/-/compare/0.2.0...0.3.0
[0.2.0]: https://gitlab.com/smspp/stochasticblock/-/compare/0.1.0...0.2.0
[0.1.0]: https://gitlab.com/smspp/stochasticblock/-/tags/0.1.0
