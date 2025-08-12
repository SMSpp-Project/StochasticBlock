/*--------------------------------------------------------------------------*/
/*-------------------- File DiscreteScenarioSet.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Class DiscreteScenarioSet that is an implementation of ScenarioGenerator
 * suited to the case where the input distribution is contained in a netCDF file
 * as a collection of vectors. Associated with the header file
 * DiscreteScenarioSet.h
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Benoit Tran \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni and Benoit Tran
 */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DiscreteScenarioSet.h"

#include "CapacitatedFacilityLocationBlock.h"
#include "Solver.h"

#include <chrono>   // For timestamp generation
#include <cstdio>   // For std::remove
#include <numeric>  // For std::accumulate (not in SMSTypedefs)

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------------- FACTORY MANAGEMENT ----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( DiscreteScenarioSet );

/*--------------------------------------------------------------------------*/
/*------------------- HELPER METHODS FOR SCENARIO REDUCTION -----------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------- HELPER METHODS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

// Get for nbScenarios
const ScenarioGenerator::ScenarioIndex &
 DiscreteScenarioSet::get_nbScenarios() const
{
  return( nbScenarios );
}

// Get for scenarioSize
const ScenarioGenerator::ScenarioSize &
 DiscreteScenarioSet::get_scenarioSize() const
{
  return( scenarioSize );
}

void DiscreteScenarioSet::empty_pool()
{
  // Clear the indexes and free memory
  scenarioIndexes.clear();
  scenarioIndexes.shrink_to_fit();
  
  // Reset pool state
  currentScenarioIndex = 0;
  poolSize = 0;
  sumPoolWeights = 0.0;
  is_initialized = false;
}

void DiscreteScenarioSet::set_poolSize( ScenarioIndex size )
{
  validate_poolSize(size);
  poolSize = size;
}

/*--------------------------------------------------------------------------*/
/*------ SCENARIO REDUCTION CONFIGURATION METHODS --------------------------*/
/*--------------------------------------------------------------------------*/

// Set scenario reduction configuration using mixed ownership model
void DiscreteScenarioSet::set_config(BlockConfig* block_config, BlockSolverConfig* solver_config)
{
  // Clean up existing configurations
  if (f_block_config) {
    delete f_block_config;
  }
  // f_solver_config will be cleaned up automatically by unique_ptr
  
  // Set the new configurations following SMS++ ownership patterns:
  // - Clone the BlockConfig (SMS++ pattern: blocks clone their configs)
  // - Take ownership of BlockSolverConfig via unique_ptr
  f_block_config = block_config ? block_config->clone() : nullptr;
  f_solver_config.reset(solver_config);
  
  // Extract and set poolSize from the BlockConfig if available
  if (block_config && block_config->f_extra_Configuration) {
    auto* poolSize_config = dynamic_cast<SimpleConfiguration<int>*>(block_config->f_extra_Configuration);
    if (poolSize_config && poolSize_config->f_value > 0) {
      poolSize = poolSize_config->f_value;
    }
  }
}

// Set scenario reduction configuration with poolSize parameter
void DiscreteScenarioSet::set_config(BlockConfig* block_config, BlockSolverConfig* solver_config, ScenarioIndex poolSize)
{
  // First set poolSize with validation
  set_poolSize(poolSize);
  
  // Then set the configurations
  set_config(block_config, solver_config);
}

// Set configuration from a Configuration object
void DiscreteScenarioSet::set_config(Configuration* config)
{
  if (!config) {
    return; // Nothing to do with null config
  }
  
  // Pattern 1: SimpleConfiguration<int> - baseline method (top poolSize by weight)
  auto* simple_poolSize = dynamic_cast<SimpleConfiguration<int>*>(config);
  if (simple_poolSize) {
    ScenarioIndex pool_size = simple_poolSize->f_value;
    if (pool_size > 0) {
      // Baseline method: just set poolSize, no CFL config needed
      poolSize = pool_size;
    }
    return;
  }
  
  // Pattern 2: SimpleConfiguration<pair<int, Configuration*>> where Configuration* is BlockSolverConfig*
  auto* poolSize_solver_pair = dynamic_cast<SimpleConfiguration<std::pair<int, Configuration*>>*>(config);
  if (poolSize_solver_pair) {
    ScenarioIndex pool_size = poolSize_solver_pair->f_value.first;
    Configuration* inner_config = poolSize_solver_pair->f_value.second;
    
    if (pool_size > 0 && inner_config) {
      // Check if inner_config is actually a BlockSolverConfig*
      auto* solver_config = dynamic_cast<BlockSolverConfig*>(inner_config);
      if (solver_config) {
        // Generate simple BlockConfig with poolSize
        auto* block_cfg = generate_default_cfl_config(pool_size);
        // Use advanced scenario reduction with provided solver
        set_config(block_cfg, solver_config, pool_size);
        // Clean up the temporary BlockConfig (set_config clones it)
        delete block_cfg;
      }
    }
    return;
  }
  
  // Pattern 3: SimpleConfiguration<pair<Configuration*, Configuration*>> where first is BlockConfig*, second is BlockSolverConfig*
  auto* config_pair = dynamic_cast<SimpleConfiguration<std::pair<Configuration*, Configuration*>>*>(config);
  if (config_pair) {
    Configuration* first_config = config_pair->f_value.first;
    Configuration* second_config = config_pair->f_value.second;
    
    if (first_config && second_config) {
      // Check if first is BlockConfig* and second is BlockSolverConfig*
      auto* block_config = dynamic_cast<BlockConfig*>(first_config);
      auto* solver_config = dynamic_cast<BlockSolverConfig*>(second_config);
      
      if (block_config && solver_config) {
        // Extract poolSize from the provided BlockConfig and use full advanced configuration
        ScenarioIndex pool_size = get_poolSize_parameter(block_config);
        // Pass the original BlockConfig (set_config will clone it)
        set_config(block_config, solver_config, pool_size);
      }
    }
    return;
  }
  
  
  // If we reach here, the configuration type is not supported
  throw std::invalid_argument(
    "Unsupported configuration type for DiscreteScenarioSet::set_config(). "
    "Supported patterns are:\n"
    "1. SimpleConfiguration<int> - baseline method (top poolSize scenarios by weight)\n"
    "2. SimpleConfiguration<pair<int, Configuration*>> - advanced with generated BlockConfig\n"
    "3. SimpleConfiguration<pair<Configuration*, Configuration*>> - full advanced configuration\n"
    "Note: The ell parameter should be set via set_ell() method or netCDF deserialization."
  );
}

// Extract poolSize parameter with validation
DiscreteScenarioSet::ScenarioIndex DiscreteScenarioSet::get_poolSize_parameter(const BlockConfig* config) const
{
  if (!config) {
    throw std::runtime_error("Invalid configuration for getting poolSize parameter");
  }
  
  // Extract poolSize from the extra configuration
  ScenarioIndex pool_size = 0;
  if (config->f_extra_Configuration) {
    auto* simple_config = dynamic_cast<SimpleConfiguration<int>*>(config->f_extra_Configuration);
    if (simple_config) {
      pool_size = simple_config->f_value;
    } else {
      throw std::runtime_error("Extra configuration is not a SimpleConfiguration<int>");
    }
  } else {
    throw std::runtime_error("No extra configuration found containing poolSize parameter");
  }
  
  // Validate the poolSize parameter
  if (pool_size == 0) {
    throw std::runtime_error("poolSize parameter must be set in configuration (value is 0)");
  }
  
  if (pool_size > nbScenarios) {
    throw std::runtime_error("Invalid poolSize parameter: " + std::to_string(pool_size) + 
                             " exceeds number of scenarios (" + std::to_string(nbScenarios) + ")");
  }
  
  return pool_size;
}

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS DiscreteScenarioSet --------------------------*/
/*--------------------------------------------------------------------------*/

void DiscreteScenarioSet::deserialize( const netCDF::NcGroup & group )
{
  // ScenarioGenerator::deserialize is pure virtual, so we start directly here
  
  // Reset state
  currentScenarioIndex = 0;
  poolSize = 0;
  sumPoolWeights = 0.0;
  is_initialized = false;

  // Clear existing data
  scenarioSet.resize(boost::extents[0][0]);
  poolProbabilities.clear();
  scenarioIndexes.clear();
  
  // Clean up any existing configuration
  if (f_block_config) {
    delete f_block_config;
    f_block_config = nullptr;
  }
  
  // f_solver_config will be reset automatically
  
  // Deserialize mandatory dimensions
  ::deserialize_dim( group , "NumberScenarios" , nbScenarios , false );
  ::deserialize_dim( group , "ScenarioSize" , scenarioSize , false );
  
  if (nbScenarios == 0) {
    throw std::invalid_argument("NumberScenarios must be positive");
  }
  if (scenarioSize == 0) {
    throw std::invalid_argument("ScenarioSize must be positive");
  }

  // Deserialize mandatory scenarios data
  scenarioSet.resize( boost::extents[ nbScenarios ][ scenarioSize ] );
  std::vector< std::size_t > sizes = { nbScenarios, scenarioSize };
  ::deserialize( group , "Scenarios" , sizes , scenarioSet , true , false );

  // Deserialize probabilities (optional, default to uniform)
  bool probsLoaded = ::deserialize(group, "poolProbabilities", nbScenarios, poolProbabilities);
  
  if (!probsLoaded) {
    // No probabilities in file, create uniform weights
    poolProbabilities.assign(nbScenarios, 1.0 / nbScenarios);
  } else if (poolProbabilities.size() != nbScenarios) {
    // Probabilities were loaded but have wrong size - this is an error
    throw std::invalid_argument("poolProbabilities size (" + std::to_string(poolProbabilities.size()) + 
                                ") does not match NumberScenarios (" + std::to_string(nbScenarios) + ")");
  }
  
  // Validate probabilities sum to approximately 1.0
  double sum = std::accumulate(poolProbabilities.begin(), poolProbabilities.end(), 0.0);
  if (std::abs(sum - 1.0) > 1e-6) {
    throw std::invalid_argument("Scenario probabilities must sum to 1.0, got: " + std::to_string(sum));
  }
  
  // Deserialize scenario reduction configuration (optional)
  try {
    netCDF::NcGroup cfgGroup = group.getGroup("ScenarioReductionConfig");
    if (!cfgGroup.isNull()) {
      // Try to read poolSize (optional)
      ScenarioIndex pool_size = 0;
      bool has_poolSize_variable = false;
      try {
        auto poolSizeVar = cfgGroup.getVar("poolSize");
        if (!poolSizeVar.isNull()) {
          poolSizeVar.getVar(&pool_size);
          has_poolSize_variable = true;
          // Note: poolSize validation will happen later during init_representative_pool
        }
      } catch (...) {
        // poolSize not found, will check BlockConfig for it
      }
      
      // Try to read ell (optional, default 2.0)
      float ell_value = 2.0f;
      try {
        auto ellVar = cfgGroup.getVar("ell");
        if (!ellVar.isNull()) {
          ellVar.getVar(&ell_value);
        }
      } catch (...) {
        // ell not found, use default
      }
      
      // Store ell as an internal variable
      this->ell = ell_value;
      
      // Try to deserialize BlockConfig
      BlockConfig* block_cfg = nullptr;
      try {
        auto blockGroup = cfgGroup.getGroup("BlockConfig");
        if (!blockGroup.isNull()) {
          auto* cfg = Configuration::new_Configuration(blockGroup);
          block_cfg = dynamic_cast<BlockConfig*>(cfg);
          if (!block_cfg) {
            delete cfg;
            throw std::runtime_error("BlockConfig deserialization failed");
          }
          // Trust the user provided appropriate config after dynamic_cast succeeded
        }
      } catch (...) {
        // No BlockConfig or invalid, will generate default if needed
        block_cfg = nullptr;
      }
      
      // Try to deserialize BlockSolverConfig
      BlockSolverConfig* solver_cfg = nullptr;
      try {
        auto solverGroup = cfgGroup.getGroup("BlockSolverConfig");
        if (!solverGroup.isNull()) {
          auto* cfg = Configuration::new_Configuration(solverGroup);
          solver_cfg = dynamic_cast<BlockSolverConfig*>(cfg);
          if (!solver_cfg) {
            delete cfg;
            throw std::runtime_error("BlockSolverConfig deserialization failed");
          }
        }
      } catch (...) {
        // No BlockSolverConfig, will generate default if needed
        solver_cfg = nullptr;
      }
      
      // Handle BlockConfig and poolSize interaction
      if (block_cfg) {
        // Check if BlockConfig has extra_Configuration as SimpleConfiguration<int>
        if (block_cfg->f_extra_Configuration) {
          auto* poolSize_config = dynamic_cast<SimpleConfiguration<int>*>(block_cfg->f_extra_Configuration);
          if (poolSize_config) {
            if (has_poolSize_variable) {
              // poolSize variable was provided, override the SimpleConfiguration value
              poolSize_config->f_value = pool_size;
            } else {
              // No poolSize variable, use SimpleConfiguration<int> value as poolSize
              pool_size = poolSize_config->f_value;
              // Note: poolSize validation will happen later during init_representative_pool
            }
          } else if (has_poolSize_variable) {
            // Extra config exists but is not SimpleConfiguration<int>, and we have poolSize
            delete block_cfg->f_extra_Configuration;
            block_cfg->f_extra_Configuration = new SimpleConfiguration<int>(pool_size);
          } else {
            // No poolSize variable and extra_config is not SimpleConfiguration<int>
            throw std::runtime_error("poolSize not found: neither as variable nor in BlockConfig's extra_Configuration");
          }
        } else if (has_poolSize_variable) {
          // No extra_Configuration but poolSize was provided, add it
          block_cfg->f_extra_Configuration = new SimpleConfiguration<int>(pool_size);
        } else {
          // No poolSize variable and no extra_Configuration
          throw std::runtime_error("poolSize not found: neither as variable nor in BlockConfig");
        }
      } else if (has_poolSize_variable) {
        // No BlockConfig provided but poolSize was given, create default one with poolSize
        block_cfg = generate_default_cfl_config(pool_size);
        if (!block_cfg) {
          throw std::runtime_error("Failed to generate default BlockConfig");
        }
      } else {
        // No BlockConfig and no poolSize variable
        throw std::runtime_error("poolSize not found: must be provided either as variable or in BlockConfig");
      }
      
      // If no solver config provided, create default
      if (!solver_cfg) {
        solver_cfg = generate_default_solver_config();
        if (!solver_cfg) {
          throw std::runtime_error("Failed to generate default BlockSolverConfig");
        }
      }
      
      // Store the configuration
      try {
        set_config(block_cfg, solver_cfg);
        poolSize = pool_size;
      } catch (const std::exception& e) {
        // Clean up allocated configs on error
        if (block_cfg) delete block_cfg;
        if (solver_cfg) delete solver_cfg;
        throw std::runtime_error(std::string("Failed to set scenario reduction config: ") + e.what());
      }
    }
  } catch (const std::exception& e) {
    // Log warning but continue - scenario reduction config is optional
    std::cerr << "Warning: Failed to load ScenarioReductionConfig: " << e.what() << std::endl;
  }
}

void DiscreteScenarioSet::serialize(netCDF::NcGroup& group) const
{
  // Always call base class serialize first
  ScenarioGenerator::serialize(group);
  
  // Serialize mandatory dimensions
  auto nbScenariosDim = group.addDim("NumberScenarios", nbScenarios);
  auto scenarioSizeDim = group.addDim("ScenarioSize", scenarioSize);
  
  // Serialize mandatory scenarios data
  auto scenariosVar = group.addVar("Scenarios", netCDF::NcDouble(), {nbScenariosDim, scenarioSizeDim});
  scenariosVar.putVar(scenarioSet.data());
  
  // Serialize probabilities
  auto probsVar = group.addVar("poolProbabilities", netCDF::NcDouble(), nbScenariosDim);
  probsVar.putVar(poolProbabilities.data());
  
  // Serialize scenario reduction configuration if it exists
  if (f_block_config && f_solver_config) {
    netCDF::NcGroup cfgGroup = group.addGroup("ScenarioReductionConfig");
    
    // Serialize poolSize if available
    if (poolSize > 0) {
      auto poolSizeVar = cfgGroup.addVar("poolSize", netCDF::ncInt);
      poolSizeVar.putVar(&poolSize);
    }
    
    // Serialize ell from internal variable
    auto ellVar = cfgGroup.addVar("ell", netCDF::ncFloat);
    ellVar.putVar(&ell);
    
    // Serialize BlockConfig
    auto blockGroup = cfgGroup.addGroup("BlockConfig");
    f_block_config->serialize(blockGroup);
    
    // Serialize BlockSolverConfig  
    auto solverGroup = cfgGroup.addGroup("BlockSolverConfig");
    f_solver_config->serialize(solverGroup);
  }
}

// Implementation for setting the seed of the pseudo-random number generator
void DiscreteScenarioSet::set_seed( unsigned long seed ) { rng.seed( seed ); }

// Draw k elements among n
/* The function generateRandomSubset draws k elements among n by use of
 * the std::shuffle function and the internal rng. The chosen indexes are
 * moved into the input variable ind. */
static void generateRandomSubset( size_t n , size_t k ,
                                  std::vector< ScenarioGenerator::ScenarioIndex > & ind ,
                                  std::mt19937 & rng )
{
  if( k > n )
    throw( std::invalid_argument( "k must be less or equal than n." ) );

  // Clear the output container
  ind.clear();
  
  // Special case: if k is 0, just return empty vector
  if (k == 0) {
    return;
  }

  // Generate ordered indexes using C++20 ranges and views
  ind.resize(k);
  
  // Create a vector with indexes 0 to n-1
  auto indexes = [n]() {
    std::vector<ScenarioGenerator::ScenarioIndex> result(n);
    std::iota(result.begin(), result.end(), 0);
    return result;
  }();
  
  // Use standard library algorithm std::sample to randomly select k items from indexes
  std::sample(indexes.begin(), indexes.end(), ind.begin(), k, rng);
}

// Initialize a pool with randomly selected scenarios
void DiscreteScenarioSet::init_random_pool(ScenarioIndex pool_size)
{
  // Clean up any existing pool
  empty_pool();
  
  // Initialize the pool parameters
  sumPoolWeights = 0.0;
  set_poolSize(pool_size);
  currentScenarioIndex = 0;
  
  // Generate random indices for the pool
  if (pool_size > 0) {
    generateRandomSubset(nbScenarios, pool_size, scenarioIndexes, rng);
    
    // Calculate the sum of weights for the selected scenarios
    sumPoolWeights = std::accumulate(scenarioIndexes.begin(), scenarioIndexes.end(), 0.0,
      [this](double sum, ScenarioIndex index) -> double {
        return sum + (index < poolProbabilities.size() ? poolProbabilities[index] : 0.0);
      });
      
    // Mark the pool as initialized
    is_initialized = true;
  }
}

// Initialize a representative pool using scenario reduction
void DiscreteScenarioSet::init_representative_pool( ScenarioIndex poolSize )
{
  // Step 1: Validate the poolSize parameter
  validate_poolSize(poolSize);
  
  // Clean up any existing pool
  empty_pool();
  
  // Check if we have a BlockSolverConfig (BlockConfig alone is not useful without solver)
  if (f_solver_config) {
    // Step 2: Create a CapacitatedFacilityLocationBlock for scenario selection
    auto cflBlock = std::make_unique<CapacitatedFacilityLocationBlock>();
    
    // Set up the scenario selection problem parameters
    ScenarioIndex n_scenarios = nbScenarios;
    ScenarioSize scenario_size = scenarioSize;
    
    // Create CFL problem data
    auto [capacities, fixed_costs, demands] = create_cfl_problem_data(n_scenarios);
    
    // Compute the transport cost matrix with ell-powered distances
    auto transport_costs = compute_transport_cost_matrix(n_scenarios, scenario_size, this->ell);
    
    // Load the CFL problem into the block
    cflBlock->load(
        n_scenarios,       // Number of facilities
        n_scenarios,       // Number of customers
        std::move(capacities),
        std::move(fixed_costs),
        std::move(demands),
        std::move(transport_costs),
        false,             // Not a balanced problem
        poolSize           // Maximum number of facilities to open
    );
    
    // Step 3.1: Apply BlockConfig if it exists, else generate minimal default
    if (f_block_config) {
      // Trust that the user provided appropriate BlockConfig for CFL
      f_block_config->apply(cflBlock.get());
    } else {
      // Generate minimal default config (mostly empty, just poolSize in extra_Configuration)
      auto* default_config = generate_default_cfl_config(poolSize);
      default_config->apply(cflBlock.get());
      delete default_config;
    }
    
    // Step 3.2: Generate abstract variables and constraints
    // With wc=7 as default, all necessary constraints are generated automatically
    #ifndef NDEBUG
    std::cout << "DEBUG [init_representative_pool]: Generating abstract variables" << std::endl;
    #endif
    cflBlock->generate_abstract_variables();
    
    #ifndef NDEBUG
    std::cout << "DEBUG [init_representative_pool]: Generating abstract constraints (using default wc=7)" << std::endl;
    std::cout << "DEBUG [init_representative_pool]: f_max_facilities = " << cflBlock->get_NMaxFacilities() << std::endl;
    #endif
    cflBlock->generate_abstract_constraints();  // Uses default wc=7
    
    // Apply BlockSolverConfig and transfer ownership to solver
    f_solver_config->apply(cflBlock.get());
    f_solver_config.release();  // Transfer ownership - we no longer own it
    
    // Get the registered solver and solve
    if (cflBlock->get_registered_solvers().empty()) {
      throw std::runtime_error("No solver registered after BlockSolverConfig::apply");
    }
    
    auto* solver = cflBlock->get_registered_solvers().front();
    if (!solver) {
      throw std::runtime_error("Failed to get solver from block");
    }
    
    #ifndef NDEBUG
    std::cout << "DEBUG [init_representative_pool]: Using solver: " << solver->classname() << std::endl;
    #endif
    
    // Solve the problem
    int status = solver->compute();
    if (status != Solver::kOK) {
      throw std::runtime_error("Solver failed with status: " + std::to_string(status));
    }
    
    // Ensure solver solution is written to block variables
    solver->get_var_solution();
    
    // Extract the solution from the CFL block
    get_selected_scenarios_from_block(cflBlock.get(), n_scenarios);
  } else {
    // No solver config - use baseline method (select top poolSize by weight)
    // BlockConfig alone is not useful for baseline selection
    apply_baseline_selection(poolSize);
  }
  
  // Update pool weights and finalize
  update_pool_weights();
  
  // Update pool size based on the selected scenarios
  poolSize = static_cast<ScenarioIndex>(scenarioIndexes.size());
  currentScenarioIndex = 0;
  
  // Mark the pool as initialized if we have scenarios
  if (!scenarioIndexes.empty()) {
    is_initialized = true;
  }
}

// Check if scenario reduction should be used
bool DiscreteScenarioSet::should_use_scenario_reduction(ScenarioIndex size) const
{
  // Skip scenario reduction for zero or single-scenario size
  if (size <= 1) return false;
  
  // Check if we have valid configurations
  if (!f_block_config || !f_solver_config) {
    return false;
  }
  
  try {
    // Try to get the poolSize parameter and validate it
    DiscreteScenarioSet::ScenarioIndex pool_size = get_poolSize_parameter(f_block_config);
    
    // If poolSize is valid, scenario reduction can be used
    return (pool_size > 0 && pool_size <= nbScenarios);
  } catch (const std::exception& e) {
    // If parameter extraction fails, don't use scenario reduction
    return false;
  }
}

// Apply scenario reduction
void DiscreteScenarioSet::apply_scenario_reduction()
{
  if (!f_block_config || !f_solver_config) {
    throw std::runtime_error("Missing scenario reduction configuration");
  }
  
  // Extract the poolSize parameter (number of scenarios to select)
  DiscreteScenarioSet::ScenarioIndex pool_size = get_poolSize_parameter(f_block_config);
  
  // DEBUG: Log poolSize parameter
  #ifndef NDEBUG
  std::cout << "DEBUG [apply_scenario_reduction]: poolSize (max scenarios to select) = " << pool_size << std::endl;
  std::cout << "DEBUG [apply_scenario_reduction]: total scenarios available = " << nbScenarios << std::endl;
  #endif
  
  // Use the internal ell parameter for distance calculations
  // (no longer extracted from BlockConfig)
  
  // Clear existing selection
  scenarioIndexes.clear();
  sumPoolWeights = 0.0;
  
  try {
    // Create a CapacitatedFacilityLocationBlock for scenario selection
    auto cflBlock = std::make_unique<CapacitatedFacilityLocationBlock>();
    
    // Set up the scenario selection problem parameters
    ScenarioIndex n_scenarios = nbScenarios;
    ScenarioSize scenario_size = scenarioSize;
    
    // Create CFL problem data using helper function
    auto [capacities, fixed_costs, demands] = create_cfl_problem_data(n_scenarios);
    
    // Compute the transport cost matrix using helper function
    auto transport_costs = compute_transport_cost_matrix(n_scenarios, scenario_size, this->ell);
    
    // Load the CFL problem into the block FIRST
    #ifndef NDEBUG
    std::cout << "DEBUG [apply_scenario_reduction]: Creating CFL block with:" << std::endl;
    std::cout << "  - n_facilities = " << n_scenarios << std::endl;
    std::cout << "  - n_customers = " << n_scenarios << std::endl;
    std::cout << "  - balanced = false" << std::endl;
    std::cout << "  - max_open_facilities = " << pool_size << std::endl;
    #endif
    
    cflBlock->load(
        n_scenarios,       // Number of facilities
        n_scenarios,       // Number of customers
        std::move(capacities),
        std::move(fixed_costs),
        std::move(demands),
        std::move(transport_costs),
        false,             // Not a balanced problem
        pool_size          // Maximum number of facilities
    );
    
    // Apply the BlockConfig AFTER loading data
    if (f_block_config) {
      f_block_config->apply(cflBlock.get());
    }
    
    // Generate abstract variables (needed for all solvers to write solution back)
    #ifndef NDEBUG
    std::cout << "DEBUG [apply_scenario_reduction]: Generating abstract variables" << std::endl;
    #endif
    cflBlock->generate_abstract_variables();
    
    // Generate abstract constraints
    // With wc=7 as default, all necessary constraints are generated automatically
    #ifndef NDEBUG
    std::cout << "DEBUG [apply_scenario_reduction]: Generating abstract constraints (using default wc=7)" << std::endl;
    std::cout << "DEBUG [apply_scenario_reduction]: f_max_facilities = " << cflBlock->get_NMaxFacilities() << std::endl;
    #endif
    cflBlock->generate_abstract_constraints();  // Uses default wc=7
    
    // Configure the solver using helper function
    Solver* solver = create_and_configure_solver(cflBlock.get(), this->ell);
    #ifndef NDEBUG
    std::cout << "DEBUG [apply_scenario_reduction]: Using solver: " << solver->classname() << std::endl;
    #endif
    
    // Solve the scenario reduction problem
    int status = solver->compute();
    #ifndef NDEBUG
    std::cout << "DEBUG [apply_scenario_reduction]: Solver status = " << status << std::endl;
    #endif
    
    // Check if the solve was successful
    if (status == Solver::kOK) {
      // Ensure solver solution is written to block variables
      solver->get_var_solution();
      
      // Extract selected scenarios using helper function
      get_selected_scenarios_from_block(cflBlock.get(), n_scenarios);
      
      // Update pool weights using helper function
      update_pool_weights();
    } else {
      throw std::runtime_error("Scenario reduction solver failed with status: " + std::to_string(status));
    }
    
  } catch (const std::exception& e) {
    // If scenario reduction fails, throw a runtime error
    throw std::runtime_error("Scenario reduction failed: " + std::string(e.what()));
  }
}


[[nodiscard]] ScenarioGenerator::Scenario DiscreteScenarioSet::get_current_scenario( void ) const
{
  if( currentScenarioIndex >= poolSize )
  {
    throw( std::out_of_range( "Current scenario index is out of range." ) );
  }

  // Make sure scenarioIndexes has the expected size
  if (scenarioIndexes.size() <= currentScenarioIndex)
  {
    throw( std::out_of_range( "scenarioIndexes is too small" ) );
  }
  
  // Make sure the index is valid
  const auto index = scenarioIndexes[currentScenarioIndex];
  if (index >= nbScenarios)
  {
    throw( std::out_of_range( "Scenario index is out of range" ) );
  }

  // Transform the scenarioIndexes[currentScenarioIndex]-th row of
  // scenarioSet into a span<const double>
  return Scenario(&scenarioSet[index][0], get_scenario_size());
}

[[nodiscard]] double DiscreteScenarioSet::get_current_scenario_probability( void ) const
{
  if( currentScenarioIndex >= poolSize )
  {
    throw( std::out_of_range( "Current scenario index is out of range." ) );
  }

  // Make sure scenarioIndexes has the expected size
  if (scenarioIndexes.size() <= currentScenarioIndex)
  {
    throw( std::out_of_range( "scenarioIndexes is too small for probability" ) );
  }
  
  // Make sure the index is valid
  const auto idx = scenarioIndexes[currentScenarioIndex];
  if (idx >= poolProbabilities.size())
  {
    throw( std::out_of_range( "Probability index is out of range" ) );
  }
  
  return (sumPoolWeights > 0.0) ? poolProbabilities[idx] / sumPoolWeights : 0.0;
}

// Implementation of the new structured binding method
[[nodiscard]] DiscreteScenarioSet::ScenarioWithProbability 
DiscreteScenarioSet::get_current_scenario_with_prob() const
{
  // Get both the scenario and probability in one call
  return {get_current_scenario(), get_current_scenario_probability()};
}

// Implementation of try_get_scenario
[[nodiscard]] std::optional<ScenarioGenerator::Scenario> 
DiscreteScenarioSet::try_get_scenario(ScenarioIndex index) const
{
  // Save current state
  if (index >= poolSize) {
    return std::nullopt;
  }
  
  try {
    if (index < scenarioIndexes.size()) {
      const auto scenarioIndex = scenarioIndexes[index];
      if (scenarioIndex < nbScenarios) {
        return Scenario(&scenarioSet[scenarioIndex][0], scenarioSize);
      }
    }
  } catch (...) {
    return std::nullopt;
  }
  
  return std::nullopt;
}

// Implementation of is_pool_initialized
[[nodiscard]] bool DiscreteScenarioSet::is_pool_initialized() const
{
  return is_initialized;
}

[[nodiscard]] bool DiscreteScenarioSet::next_scenario( void )
{
  // If poolSize is 0 or no pool is initialized, there are no scenarios to move to
  if (poolSize == 0 || !is_initialized) {
    return false;
  }
  
  if (currentScenarioIndex < poolSize - 1)
  {
    // Use prefix increment for efficiency
    ++currentScenarioIndex;
    return true; // Successfully moved to the next scenario
  }
  return false; // No more scenario in scenarioPool to move to
}

/// Implementation for retrieving the size of a scenario
[[nodiscard]] ScenarioGenerator::ScenarioSize DiscreteScenarioSet::get_scenario_size( void ) const
{
  return scenarioSize;
}

/// Concrete implementation of ScenarioGenerator
DiscreteScenarioSet::DiscreteScenarioSet() { set_seed(DEFAULT_SEED); }

/// Destructor
DiscreteScenarioSet::~DiscreteScenarioSet() {
  // Reset scenario set
  scenarioSet.resize(boost::extents[0][0]);
  
  // Clear vectors with shrink_to_fit to release memory back to the system
  scenarioIndexes.clear();
  scenarioIndexes.shrink_to_fit();
  
  poolProbabilities.clear();
  poolProbabilities.shrink_to_fit();
  
  // Release the configuration objects following SMS++ ownership patterns
  
  // BlockConfig: We always own our clone, safe to delete
  if (f_block_config) {
    delete f_block_config;
    f_block_config = nullptr;
  }
  
  // BlockSolverConfig: unique_ptr handles cleanup automatically
  // If ownership was transferred via release(), this will be nullptr and safe
}

/*--------------------------------------------------------------------------*/
/*-------------------- HELPER METHODS IMPLEMENTATION -----------------------*/
/*--------------------------------------------------------------------------*/


// Validate the poolSize parameter for scenario reduction
void DiscreteScenarioSet::validate_poolSize(ScenarioIndex poolSize) const
{
  if (poolSize == 0 || poolSize > nbScenarios) {
    throw std::invalid_argument("Invalid poolSize parameter: must be between 1 and " + 
                                std::to_string(nbScenarios));
  }
}

// Ensure scenario reduction configuration exists
void DiscreteScenarioSet::ensure_configuration_exists(ScenarioIndex poolSize)
{
  // If configuration is not initialized, create it with default values
  if (!f_block_config || !f_solver_config) {
    // Create default configurations
    auto* block_cfg = generate_default_cfl_config(poolSize);
    auto* solver_cfg = generate_default_solver_config("Dupacova");
    
    // Set the configuration
    set_config(block_cfg, solver_cfg);
    
    // Verify configuration was created
    if (!f_block_config || !f_solver_config) {
      throw std::runtime_error("Failed to create scenario reduction configuration");
    }
  } else {
    // Configuration exists, but verify it has the correct poolSize parameter
    if (f_block_config && f_block_config->f_extra_Configuration) {
      auto* poolSize_config = dynamic_cast<SimpleConfiguration<int>*>(f_block_config->f_extra_Configuration);
      if (!poolSize_config) {
        // Extra configuration exists but is not SimpleConfiguration<int>, fix it
        delete f_block_config->f_extra_Configuration;
        f_block_config->f_extra_Configuration = new SimpleConfiguration<int>(poolSize);
      }
    } else if (f_block_config) {
      // No extra configuration, add it
      f_block_config->f_extra_Configuration = new SimpleConfiguration<int>(poolSize);
    }
  }
}

// Create CFL problem data structures
std::tuple<CapacitatedFacilityLocationBlock::DVector,
           CapacitatedFacilityLocationBlock::CVector,
           CapacitatedFacilityLocationBlock::DVector>
DiscreteScenarioSet::create_cfl_problem_data(ScenarioIndex n_scenarios) const
{
  CapacitatedFacilityLocationBlock::DVector capacities(n_scenarios);
  CapacitatedFacilityLocationBlock::CVector fixed_costs(n_scenarios);
  CapacitatedFacilityLocationBlock::DVector demands(n_scenarios);
  
  // Set up the basic CFL parameters
  for (ScenarioIndex i = 0; i < n_scenarios; ++i) {
    capacities[i] = 1.0;      // Capacities = weights of the reduced distribution
    fixed_costs[i] = 0.0;     // No fixed cost in scenario reduction, but we do have a constraint on the maximal number of facilities that can (should) be opened.
    demands[i] = poolProbabilities[i];  // Demand equals probability weight
  }
  
  return std::make_tuple(std::move(capacities), std::move(fixed_costs), std::move(demands));
}

// Compute the transport cost matrix between scenarios
CapacitatedFacilityLocationBlock::CMatrix
DiscreteScenarioSet::compute_transport_cost_matrix(ScenarioIndex n_scenarios,
                                                   ScenarioSize scenario_size,
                                                   float ell) const
{
  CapacitatedFacilityLocationBlock::CMatrix transport_costs(boost::extents[n_scenarios][n_scenarios]);
  
  // Map scenarios to Eigen matrix for easier distance calculations
  Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> 
      all_scenarios(scenarioSet.data(), n_scenarios, scenario_size);
  
  // Compute the distance/transportation cost matrix
  for (ScenarioIndex i = 0; i < n_scenarios; ++i) {
    for (ScenarioIndex j = 0; j < n_scenarios; ++j) {
      if (i == j) {
        transport_costs[i][j] = 0.0;
        continue;
      }
      
      // Calculate distance between scenarios
      Eigen::VectorXd scenario_i = all_scenarios.row(i);
      Eigen::VectorXd scenario_j = all_scenarios.row(j);
      transport_costs[i][j] = compute_scenario_distance(scenario_i, scenario_j, ell);
    }
  }
  
  return transport_costs;
}

// Compute distance between two scenarios
double DiscreteScenarioSet::compute_scenario_distance(const Eigen::VectorXd& scenario1,
                                                     const Eigen::VectorXd& scenario2,
                                                     float ell) const
{
  // Transportation cost = ell-power of the euclidean norm
  return std::pow((scenario1 - scenario2).norm(), ell);
}

// Create and configure the scenario reduction solver
Solver*
DiscreteScenarioSet::create_and_configure_solver(CapacitatedFacilityLocationBlock* cflBlock,
                                                 float ell) const
{
  // Apply the BlockSolverConfig to register and configure the solver
  if (f_solver_config) {
    f_solver_config->apply(cflBlock);
    f_solver_config.release();  // Transfer ownership to solver
  }
  
  // Get the registered solver
  if (cflBlock->get_registered_solvers().empty()) {
    throw std::runtime_error("No solver registered to the block after BlockSolverConfig::apply");
  }
  
  auto* base_solver = cflBlock->get_registered_solvers().front();
  if (!base_solver) {
    throw std::runtime_error("Failed to get solver from block");
  }
  
  // Return the solver pointer (block owns it)
  return base_solver;
}

// Extract selected scenarios from solver results
void DiscreteScenarioSet::get_selected_scenarios_from_block(const CapacitatedFacilityLocationBlock* cflBlock,
                                                          ScenarioIndex n_scenarios)
{
  #ifndef NDEBUG
  std::cout << "DEBUG [get_selected_scenarios_from_block]: Starting extraction from CFL block" << std::endl;
  #endif
  
  // Clear existing selection
  scenarioIndexes.clear();
    
  // Read variable values directly from the block variables
  for (ScenarioIndex i = 0; i < n_scenarios; ++i) {
    // Get the y variable for facility i
    const auto* y_var = cflBlock->get_y(i);
    if (!y_var) {
      throw std::runtime_error("Failed to get y variable for facility " + std::to_string(i));
    }
    
    // Get value from the variable directly
    double y_value = y_var->get_value();
    #ifndef NDEBUG
    std::cout << "DEBUG [get_selected_scenarios_from_block]: y[" << i << "] = " << y_value << std::endl;
    #endif
    
    // Check if facility i is open (y[i] > 0.5)
    if (y_value > 0.5) {
      scenarioIndexes.push_back(i);
    }
  }
  
  #ifndef NDEBUG
  std::cout << "DEBUG [get_selected_scenarios_from_block]: Selected " << scenarioIndexes.size() 
            << " scenarios: ";
  for (auto idx : scenarioIndexes) std::cout << idx << " ";
  std::cout << std::endl;
  #endif
  
  // If no scenarios were selected, throw an error
  if (scenarioIndexes.empty()) {
    throw std::runtime_error("No scenarios selected by the reduction algorithm");
  }
}

// Apply baseline selection method
void DiscreteScenarioSet::apply_baseline_selection(ScenarioIndex poolSize)
{
  // Create pairs of (index, weight) for sorting
  std::vector<std::pair<ScenarioIndex, double>> indexed_weights;
  indexed_weights.reserve(nbScenarios);
  
  for (ScenarioIndex i = 0; i < nbScenarios; ++i) {
    double weight = (i < poolProbabilities.size()) ? poolProbabilities[i] : 1.0 / nbScenarios;
    indexed_weights.emplace_back(i, weight);
  }
  
  // Sort by weight in descending order (highest weights first)
  std::sort(indexed_weights.begin(), indexed_weights.end(),
            [](const auto& a, const auto& b) { return a.second > b.second; });
  
  // Select the top poolSize scenarios
  scenarioIndexes.clear();
  scenarioIndexes.reserve(poolSize);
  
  for (ScenarioIndex i = 0; i < poolSize && i < nbScenarios; ++i) {
    scenarioIndexes.push_back(indexed_weights[i].first);
  }
  
  // Sort the selected indices for consistency
  std::sort(scenarioIndexes.begin(), scenarioIndexes.end());
  
  #ifndef NDEBUG
  std::cout << "DEBUG [apply_baseline_selection]: Selected " << scenarioIndexes.size() 
            << " scenarios using baseline method (top weights)" << std::endl;
  #endif
}

// Update pool weights after scenario selection
void DiscreteScenarioSet::update_pool_weights()
{
  // Calculate sum of weights of selected scenarios for normalization
  sumPoolWeights = 0.0;
  for (const auto& idx : scenarioIndexes) {
    sumPoolWeights += poolProbabilities[idx];
  }
}



// Generate default BlockConfig for CFL
BlockConfig* DiscreteScenarioSet::generate_default_cfl_config(ScenarioIndex poolSize) const
{
  auto* config = new BlockConfig(false);  // not differential
  
  // Set poolSize in extra_Configuration
  config->f_extra_Configuration = new SimpleConfiguration<int>(poolSize);
  
  return config;
}

// Generate default BlockSolverConfig
BlockSolverConfig* DiscreteScenarioSet::generate_default_solver_config(const std::string& algorithm) const
{
  auto* config = new BlockSolverConfig(true);  // differential
  
  // Add a default solver configuration
  // Note: ScenarioReductionSolver is used as default, but user can provide any suitable solver
  config->add_ComputeConfig("ScenarioReductionSolver", nullptr);
  
  return config;
}

/*--------------------------------------------------------------------------*/
/*------------------ End file DiscreteScenarioSet.cpp ----------------------*/
/*--------------------------------------------------------------------------*/