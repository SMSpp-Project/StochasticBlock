# StochasticBlock

SMS++ `Block` that provides the "plumbing" for representing (possibly,
among other things) stochastic optimization models.

The `StochasticBlock` class is meant to represent a `Block` whose data
changes in a "scenario way", i.e., whole swathes of the data change all in
one blow. This is typical (among other things) in stochastic optimization
models. While `StochasticBlock` itself focuses on the mechanics of applying
scenario data to modify its inner Block, it is designed to work seamlessly
with `ScenarioGenerator` classes that handle probability distributions and
scenario generation. A `StochasticBlock` is characterized by the following:

- **Inner Block**: A single inner `Block` whose data changes; this can be
  any `:Block`

- **Scenario Data**: Some data of the inner `Block` is subject to changes,
  and the value for this data is represented by a vector of double; an
  instance of this vector is called a *scenario* for the data (a common
  term in stochastic optimization)

- **Data Mappings**: A set of `DataMapping` objects that identify the data
  in the inner  `Block` that is subject to change and specify how to modify
  this data. The inner `Block` may have different pieces of changing data,
  located in different sub-`Block`s; each `DataMapping` represents one of
  these pieces of data

- **Scenario Application**: A `set_data()` method which takes a scenario (a
  vector) as  parameter and sets the data of the inner `Block` according to
  its set of `DataMapping`

- **ScenarioGenerator Integration**: Designed to work alongside
  `ScenarioGenerator` classes  (such as `DiscreteScenarioSet`) that manage
  probability distributions and provide scenarios. The separation of
  concerns is deliberate: `ScenarioGenerator` handles the stochastic aspects
  (probabilities, sampling, reduction), while `StochasticBlock` handles the
  application of scenario data to the optimization model

## Usage Patterns

### Typical Workflow

1. Create a `StochasticBlock` with an inner deterministic block

2. Define data mappings to specify which data varies per scenario

3. Create a `ScenarioGenerator` (e.g., `DiscreteScenarioSet`) to provide
   scenario data

4. Iterate through scenarios from the generator, applying each via
   `set_data()`

5. Solve or process the block for each scenario realization

## Creating a StochasticBlock

### Parametric Construction

You can create a `StochasticBlock` programmatically:

```cpp
// Create the inner deterministic block
auto inner_block = new MyDeterministicBlock();

// Create the StochasticBlock wrapper
auto stochastic_block = new StochasticBlock(nullptr, inner_block);

// Add data mappings to specify which data changes per scenario
// SimpleDataMappingBase is the usual class for data mappings
auto mapping1 = std::make_unique<SimpleDataMapping<double>>(
    inner_block,                      // block containing the data
    "cost_vector",                    // name/identifier of the data
    10                                // size of data
    );
stochastic_block->add_data_mapping(std::move(mapping1));

// Option 1: Apply scenarios manually
std::vector<double> scenario_data = {/* ... */};
stochastic_block->set_data(scenario_data);

// Option 2: Use a ScenarioGenerator (kept separate from StochasticBlock)
auto* scenario_gen = new DiscreteScenarioSet();
// ... configure and load scenarios into scenario_gen ...

// Iterate through scenarios from the generator
do {
    auto scenario = scenario_gen->get_current_scenario();
    stochastic_block->set_data(scenario);  // Direct Scenario object support
    // ... solve or process the block with this scenario ...
    } while (scenario_gen->next_scenario());
```

### Deserialization from netCDF

The `StochasticBlock` can be loaded from a netCDF file using the
`deserialize()` method. The expected structure is:

```
StochasticBlock Group
├── Attribute: type = "StochasticBlock"          [Required]
│
├── Group: Block/                                 [Optional]
│   ├── Attribute: type = "<BlockType>"          (e.g., "UCBlock", "MCFBlock")
│   └── ... block-specific data ...
│
├── Dimension: NumberDataMappings = N            [Optional]
│
├── Group: DataMapping_0/                        [Optional, if NumberDataMappings > 0]
│   ├── Attribute: type = "<MappingType>"        (e.g., "SimpleDataMapping")
│   ├── Attribute: offset = <int>                (position in scenario vector)
│   ├── Attribute: size = <int>                  (number of elements)
│   ├── Attribute: block_path = "<string>"       (path to sub-block)
│   └── Attribute: data_name = "<string>"        (name of data field)
│
├── Group: DataMapping_1/                        [Optional]
│   └── ... same structure as DataMapping_0 ...
│
└── Group: DataMapping_N/                        [Additional mappings as needed]
```

**Key points:**

- The inner `Block` group is optional - it can be provided separately via
  `set_inner_block()`

- Data mappings are optional but typically needed for the block to be useful

- Each data mapping tells the StochasticBlock which parts of the scenario
  vector correspond to which data in the inner block

- ScenarioGenerators (like DiscreteScenarioSet) are kept separate and used
  alongside StochasticBlock

## DiscreteScenarioSet

The `DiscreteScenarioSet` class is a concrete implementation of
`ScenarioGenerator` that manages discrete probability distributions
represented as collections of scenario vectors.

### Key Features

- **Data Management**: Scenarios are loaded from netCDF files and stored
  internally as a `boost::multi_array<double, 2>` where each row represents
  a scenario vector

- **Probability Weights**: Supports both uniform and weighted probability
  distributions over scenarios

- **Scenario Selection Methods**:

  * **Random Selection** (`init_random_pool()`): Randomly samples scenarios
    from the full set using weighted sampling

  * **Representative Selection** (`init_representative_pool()`): Selects a
    representative subset that minimizes the Wasserstein distance (when
	configured with solver) or selects highest-weight scenarios (baseline
	method)

### Usage Pattern

```cpp
// 1. Load scenarios from netCDF
DiscreteScenarioSet scenario_set;
netCDF::NcFile file("scenarios.nc", netCDF::NcFile::read);
scenario_set.deserialize(file.getGroup("Scenarios"));

// 2. Configure scenario reduction (optional)
scenario_set.set_config(block_config, solver_config);

// 3. Initialize pool
scenario_set.init_random_pool(100);        // Random sampling
// OR
scenario_set.init_representative_pool(50); // Representative selection

// 4. Iterate through scenarios
do {
    auto scenario = scenario_set.get_current_scenario();
    double probability = scenario_set.get_current_scenario_probability();
    // Process scenario...
    } while (scenario_set.next_scenario());
```

### Scenario Reduction

When configured with appropriate solvers, `DiscreteScenarioSet` can perform
optimization-based scenario reduction by formulating the problem as a
`CapacitatedFacilityLocationBlock` instance. This minimizes the Wasserstein
distance between the original and reduced distributions.

The class offers **three ways** to configure scenario reduction, allowing
you to choose the approach that best fits your workflow:

#### Configuration Methods

1. **Direct Configuration with Separate Parameters**:
   ```cpp
   auto* block_config = new BlockConfig();
   auto* solver_config = new BlockSolverConfig();
   // Configure parameters as needed...
   
   // Option A: Set config without pool size (specify later)
   scenario_set.set_config(block_config, solver_config);
   
   // Option B: Set config with pool size together
   scenario_set.set_config(block_config, solver_config, 50);
   ```

2. **Configuration via netCDF Deserialization**:

   - Configuration can be loaded automatically from netCDF files

   - The file should contain a `ScenarioReductionConfig` group with
     serialized BlockConfig and BlockSolverConfig

   - This enables saving and reusing configurations across runs

   ```cpp
   // Configuration loaded automatically during deserialization
   scenario_set.deserialize(nc_group);  // Loads data AND configuration if present
   ```

3. **Configuration Objects (SMS++ Style)**:

   ```cpp
   // Simple: Just pool size (uses baseline method)
   auto config = new SimpleConfiguration<int>(50);
   scenario_set.set_config(config);
   
   // Advanced: Pool size + Solver configuration
   auto config = new SimpleConfiguration<pair<int, Configuration*>>(
       make_pair(50, solver_config)
   );
   scenario_set.set_config(config);
   
   // Full: Complete Block and Solver configuration
   auto config = new SimpleConfiguration<pair<Configuration*, Configuration*>>(
       make_pair(block_config, solver_config)
   );
   scenario_set.set_config(config);
   ```

Without BlockSolverConfig-uration, `init_representative_pool()` falls back
to selecting scenarios with the highest probability weights.


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

- **Benoît Tran**  
  Dipartimento di Informatica  
  Università di Pisa

### Contributors

- **Antonio Frangioni**  
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
