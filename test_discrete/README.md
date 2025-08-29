# DiscreteScenarioSet Tests

Test suite for the DiscreteScenarioSet class.

## Overview

Tests scenario management and reduction functionality including:
- Random and representative scenario pool selection
- ScenarioReductionSolver algorithms (Dupacova, BestFit, FirstFit)
- MILPSolver support (CPLEX, Gurobi, SCIP, HiGHS)
- Configuration handling and error cases

## Running Tests

```bash
# Build the test suite
make

# Run all tests
./StochasticBlock_test_discretescenarioset

# List all available tests
./StochasticBlock_test_discretescenarioset -l

# Run specific test by name
./StochasticBlock_test_discretescenarioset -t "Scenario Reduction - MILP"

# Run with verbose output
./StochasticBlock_test_discretescenarioset -v

# Show help
./StochasticBlock_test_discretescenarioset -h
```

## Available Tests

- Basic Loading
- Invalid k=0
- Invalid k>nbScenarios  
- Valid k values
- Random Pool
- Scenario Reduction - Dupacova
- Scenario Reduction - BestFit
- Scenario Reduction - MILP
- Invalid Solver Config
- Algorithm Comparison

## Requirements

- SMS++ framework
- StochasticBlock module
- CapacitatedFacilityLocationBlock module
- netCDF library
- Optional: MILPSolver with backend (CPLEX, Gurobi, SCIP, or HiGHS)