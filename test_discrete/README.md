# DiscreteScenarioSet Test Suite

Test suite for the `DiscreteScenarioSet` class, covering scenario management, reduction algorithms, and configuration patterns.

## Overview

This test suite validates the functionality of `DiscreteScenarioSet`, which manages discrete probability distributions for stochastic optimization. The tests cover:

- **Scenario Management**: Loading, storing, and accessing scenario data
- **Pool Selection**: Random sampling and representative selection methods
- **Scenario Reduction**: Optimization-based algorithms to find representative subsets
- **Configuration**: Multiple configuration patterns and serialization
- **Error Handling**: Edge cases and invalid parameter handling

## Building the Tests

### Using Make
```bash
cd test_discrete
make
```

### Using CMake
```bash
cd test_discrete
mkdir build && cd build
cmake ..
make
```

## Running the Tests

### Basic Usage
```bash
# Run all tests
./test_discrete

# Run with verbose output (shows timing and details)
./test_discrete -v

# Run with minimal output
./test_discrete -q
```

### Command Line Options
- `-v, --verbose`: Enable verbose output with timing information
- `-q, --quiet`: Minimal output (only failures)
- `-h, --help`: Display help message

## Requirements

### Required Dependencies
- SMS++ core library
- StochasticBlock
- CapacitatedFacilityLocationBlock
- netCDF-C++ library
- Eigen3
- Boost (multi_array)

### Optional Dependencies
For full scenario reduction testing:
- MILPSolver module with at least one backend:
  - CPLEX (commercial)
  - Gurobi (commercial)
  - SCIP (academic)
  - HiGHS (open-source)

### Debug Build
```bash
# Enable debug symbols and assertions
make clean
make SMS_DEBUG=1
```

## Contributing

When adding new tests:
1. Follow the existing test structure (Test N format)
2. Clean up any temporary files created
3. Update this README

## License

See the main StochasticBlock [LICENSE](../LICENSE) file.