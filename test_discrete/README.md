# DiscreteScenarioSet Tests

Tests for the `DiscreteScenarioSet` class, part of the StochasticBlock module.

## Overview

The `DiscreteScenarioSet` provides two pool types with multiple selection methods:

1. **Discrete Pool**:
   - Random selection: Uniform random subset of original scenarios
   - Scenario reduction: Optimal subset selection using Wasserstein distance minimization (requires CapacitatedFacilityLocationBlock)

2. **Continuous Pool**: 
   - Representative scenarios via k-means clustering

The implementation uses modern C++ features including range-based algorithms, structured bindings, std::optional, constexpr, lambda expressions, [[nodiscard]] attributes, and std::variant for pool types. The scenario reduction functionality integrates with SMS++ Configuration system for flexible algorithm customization.

## Running Tests

```bash
# Build the test suite
make

# Run all tests
./test_discretescenarioset

# Run specific tests (by number)
./test_discretescenarioset 1 2 3

# Run with verbose output
./test_discretescenarioset -v

# Clean up temporary files without running tests
./test_discretescenarioset --clean
```

## Available Tests

1. Basic Deserialization
2. Input Validation
3. Scenario Access
4. Random Seed
5. Empty State
6. Pool Switching
7. Probability Distribution
8. Edge Cases
9. Memory Management
10. Large Scenario Set (Scalability)
11. Continuous Pool
12. Configuration Integration
13. Scenario Reduction Functionality

Note: Tests 10 and 13 can take more time to run due to large dataset processing and optimization algorithms.