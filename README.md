# StochasticBlock

SMS++ `Block` that provides the "plumbing" for representing (possibly,
among other things) stochastic optimization models.

The `StochasticBlock` class is meant to represent a `Block` whose data
changes in a "scenario way", i.e., whole swathes of the data change all in
one blow. This is typical (among other things) in stochastic optimization
models, although `StochasticBlock` per se does not have any direct notion
of stochasticity. Indeed, a `StochasticBlock` is characterized by the
following:

- a single inner `Block` whose data changes; this can be any `:Block`

- some data of the inner `Block` is subject to changes, and the value
  for this data is represented by a vector of double; a an instance of
  this vector is called a *scenario* for the data (a common term in
  stochastic optimization), but note that this scenario information is
  not stored in the `StochasticBlock`

- `StochasticBlock` has a set of `DataMapping`. This is used to both
  identify the data in the inner `Block` that is subject to change and
  to modify the values of this data. The inner `Block` may have different
  pieces of changing data, located in different sub-`Block` of its; a
  single `DataMapping` is meant to represent one of these pieces data.

- `StochasticBlock` has a `set_data()` method which takes has a scenario
  (a vector) as parameter and sets the data of the inner `Block` according
  to its set of `DataMapping`

### DiscreteScenarioSet

The `DiscreteScenarioSet` class extends `ScenarioGenerator` to manage discrete sets of scenarios loaded from data files or generated programmatically. Key features include:

- **Scenario Pool Management**: Supports both random sampling and representative pool selection
- **Serialization**: Read/write scenarios from/to netCDF files with associated probabilities
- **Scenario Reduction**: When built with `CapacitatedFacilityLocationBlock`, provides advanced scenario reduction capabilities (see below)

### Scenario Reduction

The `DiscreteScenarioSet` class provides scenario reduction capabilities when built with the `CapacitatedFacilityLocationBlock` dependency. This feature allows for reducing a large set of scenarios to a smaller representative subset while minimizing the Wasserstein distance between the original and reduced distributions.

#### Available Algorithms

- **Dupacova algorithm**: The default scenario reduction method that uses a forward selection approach to minimize Wasserstein distance
- **BestFit**: An alternative algorithm that uses local search with best improvement selection at each step
- **FirstFit**: A simpler, faster algorithm that selects the first satisfactory improvement

#### Configuration

Scenario reduction can be configured in three ways:

1. **Manual Configuration**: Create BlockConfig and BlockSolverConfig objects:
   ```cpp
   // Create BlockConfig with k parameter
   auto* block_config = new BlockConfig();
   block_config->f_extra_Configuration = new SimpleConfiguration<int>(10);  // k=10
   block_config->f_static_variables_Configuration = new SimpleConfiguration<double>(2.0);  // ell=2.0
   
   // Create BlockSolverConfig for solver selection
   auto* solver_config = new BlockSolverConfig(true);  // differential mode
   solver_config->add_ComputeConfig("ScenarioReductionSolver", nullptr);  // For Dupacova algorithm
   // Or for MILP solver:
   // solver_config->add_ComputeConfig("CPXMILPSolver", compute_config);
   
   // Apply to scenario set
   discreteScenarioSet.set_scenario_reduction_config(block_config, solver_config);
   ```

2. **From netCDF File**: The DiscreteScenarioSet can deserialize previously saved configurations from netCDF files. When a `ScenarioReductionConfig` group is present in the file, it will be loaded automatically. Note that this requires the configuration to have been previously serialized using DiscreteScenarioSet's serialize method - you cannot manually create this structure in the netCDF file as it contains serialized BlockConfig and BlockSolverConfig objects. Will be fixed soonish.

3. **Alternative Configuration Approach**: Using nested configuration objects:
   ```cpp
   auto config = new BlockConfig(false);  // Not differential
   auto cflConfig = new BlockConfig(false);
   cflConfig->add_conf_prop("k", 5);
   cflConfig->add_conf_prop("ell", 2.0f);
   config->add_R3_BlockConfig("CFLConfig", cflConfig);
   
   auto solverConfig = new BlockSolverConfig(false);  // Not differential
   solverConfig->add_conf_prop("algorithm", std::string("Dupacova"));
   config->add_Solver_Configuration("SolverConfig", solverConfig);
   
   discreteScenarioSet.set_scenario_reduction_config(config);
   ```

When the `CapacitatedFacilityLocationBlock` is not available, the class automatically falls back to a simpler selection method based on scenario probabilities.


## Getting started

These instructions will let you build StochasticBlock on your system.

### Requirements

- [SMS++ core library](https://gitlab.com/smspp/smspp)
- [CapacitatedFacilityLocationBlock](https://gitlab.com/smspp/capacitatedfacilitylocationblock) (optional, enables scenario reduction functionality)

### Build and install with CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
cmake --build .
```

The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration).

Optionally, install the library in the system with:

```sh
cmake --install .
```

### Usage with CMake

After the library is built, you can use it in your CMake project with:

```cmake
find_package(StochasticBlock)
target_link_libraries(<my_target> SMS++::StochasticBlock)
```

### Running the tests with CMake

Some unit tests will be built with the library.
Launch `ctest` from the build directory to run them.
To disable them, set the option `BUILD_TESTING` to `OFF`.

### Build and install with makefiles

Carefully hand-crafted makefiles have also been developed for those unwilling
to use CMake. Makefiles build the executable in-source (in the same directory
tree where the code is) as opposed to out-of-source (in the copy of the
directory tree constructed in the build/ folder) and therefore it is more
convenient when having to recompile often, such as when developing/debugging
a new module, as opposed to the compile-and-forget usage envisioned by CMake.

Each executable using `StochasticBlock` has to include a "main makefile" of
the module, which typically is either [makefile-c](makefile-c) including all
necessary libraries comprised the "core SMS++" one, or
[makefile-s](makefile-s) including all necessary libraries but not the "core
SMS++" one (for the common case in which this is used together with other
modules that already include them). One relevant case is the [unit
tests](test/test.cpp). The makefiles in turn recursively include all the
required other makefiles, hence one should only need to edit the "main
makefile" for compilation type (C++ compiler and its options) and it all
should be good to go. In case some of the external libraries are not at
their default location, it should only be necessary to create the
`../extlib/makefile-paths` out of the `extlib/makefile-default-paths-*` for
your OS `*` and edit the relevant bits (commenting out all the rest).

Check the [SMS++ installation wiki](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration#location-of-required-libraries)
for further details.


## Getting help

If you need support, you want to submit bugs or propose a new feature, you can
[open a new issue](https://gitlab.com/smspp/stochasticblock/-/issues/new).


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.

## Authors

### Lead Authors

- **Rafael Durbano Lobato**  
  Dipartimento di Informatica  
  Università di Pisa

### Contributors

- **Benoît Tran**  
  Dipartimento di Informatica  
  Università di Pisa


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details.


## Disclaimer

The code is currently provided free of charge under an open-source license.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.
