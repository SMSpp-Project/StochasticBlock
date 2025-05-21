# DiscreteScenarioSet Tests

Tests for the `DiscreteScenarioSet` class, part of the StochasticBlock module.

## Overview

The `DiscreteScenarioSet` provides functionality for scenario management and reduction with multiple selection methods:

1. **Scenario Pool Management**:
   - Random selection: Uniform random subset of original scenarios
   - Scenario reduction: Optimal subset selection using Wasserstein distance minimization
   - Configuration-based initialization and customization

2. **Advanced Features**: 
   - Clustered test data generation for realistic scenario distributions
   - Wasserstein distance calculation for quantifying scenario reduction quality
   - Exception safety and robustness in edge cases

The implementation integrates with the SMS++ Configuration system for flexible algorithm customization and supports multiple scenario reduction algorithms.

## Running Tests

```bash
# Build the test suite
make

# Run all tests
./StochasticBlock_test_discretescenarioset

# Run specific tests (by number)
./StochasticBlock_test_discretescenarioset 1 8 9

# Run with verbose output
./StochasticBlock_test_discretescenarioset -v

# Run a reduced set of tests (faster)
./StochasticBlock_test_discretescenarioset -q
```

## Available Tests

1. Basic Deserialization: Tests basic file loading and scenario deserialization
2. Random Pool Initialization: Tests creating a random subset of scenarios
3. Scenario Reduction Configuration: Tests configuration management
8. Robust Configuration Handling: Tests loading and using embedded configurations
9. Algorithm Comparison with Exception Safety: Tests multiple reduction algorithms with error handling
10. Deterministic Results: Tests consistency of results with fixed seeds
11. Stress Test with Robustness: Tests performance and stability with large datasets
12. Edge Cases Handling: Tests behavior in boundary conditions

Note:
- Tests 1-3 are basic functionality tests that run quickly
- Tests 8-12 are more comprehensive tests added from the improved test suite
- Tests 10-11 may take more time to run and are excluded from the quick test mode (-q)