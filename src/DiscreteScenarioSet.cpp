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

#include <span>  // Ensure span is available before other includes
#include "DiscreteScenarioSet.h" 
#include "CapacitatedFacilityLocationBlock.h"
#include "ScenarioReductionSolver.h"
#include "Block.h" 
#include "Solver.h"
#include "BlockSolverConfig.h" // For BlockSolverConfig class

using namespace SMSpp_di_unipi_it;

// Needed for unique temporary filename generation and file operations
#include <chrono>
#include <cstdio>  // for std::remove
#include <numeric> // for std::accumulate

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
  if( size > nbScenarios ) // indirectly protects if input size is negative
    throw( std::out_of_range( "The desired sample size is greater than "
                              "the number of available number of different "
                              "scenarios." ) );
  poolSize = size;
}

/*--------------------------------------------------------------------------*/
/*------ SCENARIO REDUCTION CONFIGURATION METHODS --------------------------*/
/*--------------------------------------------------------------------------*/

// Get BlockConfig for scenario reduction
BlockConfig* DiscreteScenarioSet::get_scenario_reduction_block_config() const
{
  return f_scenario_reduction_config.first;
}

// Get solver configuration for scenario reduction
BlockSolverConfig* DiscreteScenarioSet::get_scenario_reduction_solver_config() const
{
  return f_scenario_reduction_config.second;
}

// Set scenario reduction configuration
void DiscreteScenarioSet::set_scenario_reduction_config(BlockConfig* block_config, BlockSolverConfig* solver_config)
{
  // Validate the configurations
  if (block_config && solver_config) {
    // Validate BlockSolverConfig contains a valid solver name
    const auto& solver_names = solver_config->get_SolverNames();
    if (solver_names.empty()) {
      throw std::invalid_argument("BlockSolverConfig must specify at least one solver");
    }
    
    // Check if the solver is supported
    const std::string& solver_name = solver_names[0];
    bool is_scenario_reduction_solver = (solver_name == "ScenarioReductionSolver");
    bool is_milp_solver = (solver_name == "CPXMILPSolver" || 
                          solver_name == "GRBMILPSolver" || 
                          solver_name == "SCIPMILPSolver" || 
                          solver_name == "HiGHSMILPSolver");
    
    if (!is_scenario_reduction_solver && !is_milp_solver) {
      throw std::invalid_argument("Unsupported solver for scenario reduction: " + solver_name + 
                                ". Supported solvers are: ScenarioReductionSolver, CPXMILPSolver, GRBMILPSolver, SCIPMILPSolver, HiGHSMILPSolver");
    }
    
    // For ScenarioReductionSolver, validate algorithm parameter if present
    if (is_scenario_reduction_solver) {
      const auto& solver_configs = solver_config->get_SolverConfigs();
      if (!solver_configs.empty() && solver_configs[0]) {
        // Check if algorithm parameter is valid
        auto* compute_config = solver_configs[0];
        // Find strAlgorithm in str_pars vector
        for (const auto& [param_name, param_value] : compute_config->str_pars) {
          if (param_name == "strAlgorithm") {
            const std::string& algorithm = param_value;
            if (algorithm != "Dupacova" && algorithm != "BestFit" && algorithm != "FirstFit") {
              throw std::invalid_argument("Invalid algorithm for ScenarioReductionSolver: " + algorithm + 
                                        ". Valid options are: Dupacova, BestFit, FirstFit");
            }
            break;
          }
        }
      }
    }
  }
  
  // Clean up existing configurations if they are different from the new ones
  if (f_scenario_reduction_config.first && f_scenario_reduction_config.first != block_config) {
    delete f_scenario_reduction_config.first;
  }
  
  if (f_scenario_reduction_config.second && f_scenario_reduction_config.second != solver_config) {
    delete f_scenario_reduction_config.second;
  }
  
  // Set the new configurations
  f_scenario_reduction_config.first = block_config;
  f_scenario_reduction_config.second = solver_config;
  
  // Extract and set k_value from the BlockConfig if available
  if (block_config && block_config->f_extra_Configuration) {
    auto* k_config = dynamic_cast<SimpleConfiguration<int>*>(block_config->f_extra_Configuration);
    if (k_config && k_config->f_value > 0) {
      k_value = k_config->f_value;
    }
  }
}

// Set scenario reduction configuration with k parameter
void DiscreteScenarioSet::set_scenario_reduction_config(BlockConfig* block_config, BlockSolverConfig* solver_config, ScenarioIndex k)
{
  // First set k_value with validation
  set_k_value(k);
  
  // Then set the configurations
  set_scenario_reduction_config(block_config, solver_config);
}

// Set configuration from a Configuration object
void DiscreteScenarioSet::set_config(Configuration* config)
{
  if (!config) {
    return; // Nothing to do with null config
  }
  
  // For scenario reduction, we expect the Configuration to directly contain
  // k, BlockConfig, and BlockSolverConfig, not wrapped in ComputeConfig
  
  // First check if this is a BlockConfig (contains scenario reduction settings)
  auto* block_config = dynamic_cast<BlockConfig*>(config);
  if (block_config) {
    // Extract k from the extra configuration
    ScenarioIndex k = 0;
    if (block_config->f_extra_Configuration) {
      auto* k_config = dynamic_cast<SimpleConfiguration<int>*>(block_config->f_extra_Configuration);
      if (k_config) {
        k = k_config->f_value;
      }
    }
    
    if (k > 0) {
      // Create solver config if not present (use default)
      auto* solver_cfg = generate_default_solver_config("Dupacova");
      set_scenario_reduction_config(block_config->clone(), solver_cfg, k);
      
      // If scenarios are already loaded, apply reduction immediately
      if (nbScenarios > 0) {
        init_representative_pool(k);
      } else {
        // Store k for later application
        pending_k_value = k;
      }
    }
    return;
  }
  
  // Otherwise, check if it's a generic Configuration with sub-configs
  // This would be the case when loaded from netCDF with ScenarioReductionConfig group
  
  // Try to extract k parameter directly
  ScenarioIndex k = 0;
  double ell = DEFAULT_ELL_VALUE;
  
  // Check for SimpleConfiguration<int> (just k)
  auto* k_config = dynamic_cast<SimpleConfiguration<int>*>(config);
  if (k_config) {
    k = k_config->f_value;
  } else {
    // Check for SimpleConfiguration<vector<double>> format: [k, ell, ...]
    auto* vec_config = dynamic_cast<SimpleConfiguration<std::vector<double>>*>(config);
    if (vec_config && !vec_config->f_value.empty()) {
      k = static_cast<ScenarioIndex>(vec_config->f_value[0]);
      if (vec_config->f_value.size() > 1) {
        ell = vec_config->f_value[1];
      }
    }
  }
  
  // If we found k, create default configs and apply
  if (k > 0) {
    auto* block_cfg = generate_default_cfl_config(k, ell);
    auto* solver_cfg = generate_default_solver_config("Dupacova");
    
    set_scenario_reduction_config(block_cfg, solver_cfg, k);
    
    // If scenarios are already loaded, apply reduction immediately
    if (nbScenarios > 0) {
      init_representative_pool(k);
    } else {
      // Store k for later application
      pending_k_value = k;
    }
  }
}

// Extract k parameter with validation
DiscreteScenarioSet::ScenarioIndex DiscreteScenarioSet::get_k_parameter(const BlockConfig* config) const
{
  if (!config) {
    throw std::runtime_error("Invalid configuration for getting k parameter");
  }
  
  // Extract k from the extra configuration
  ScenarioIndex k = 0;
  if (config->f_extra_Configuration) {
    auto* simple_config = dynamic_cast<SimpleConfiguration<int>*>(config->f_extra_Configuration);
    if (simple_config) {
      k = simple_config->f_value;
    } else {
      throw std::runtime_error("Extra configuration is not a SimpleConfiguration<int>");
    }
  } else {
    throw std::runtime_error("No extra configuration found containing k parameter");
  }
  
  // Validate the k parameter
  if (k == 0) {
    throw std::runtime_error("k parameter must be set in configuration (value is 0)");
  }
  
  if (k > nbScenarios) {
    throw std::runtime_error("Invalid k parameter: " + std::to_string(k) + 
                             " exceeds number of scenarios (" + std::to_string(nbScenarios) + ")");
  }
  
  return k;
}

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS DiscreteScenarioSet --------------------------*/
/*--------------------------------------------------------------------------*/

void DiscreteScenarioSet::deserialize( const netCDF::NcGroup & group )
{
  // ScenarioGenerator::deserialize is pure virtual, so we start directly here
  
  // Reset state to properly handle multiple deserializations
  currentScenarioIndex = 0;
  poolSize = 0;
  sumPoolWeights = 0.0;
  is_initialized = false;

  // Clear existing data
  scenarioSet.resize(boost::extents[0][0]);
  poolProbabilities.clear();
  scenarioIndexes.clear();
  
  // Clean up any existing configuration
  if (f_scenario_reduction_config.first) {
    delete f_scenario_reduction_config.first;
    f_scenario_reduction_config.first = nullptr;
  }
  
  if (f_scenario_reduction_config.second) {
    delete f_scenario_reduction_config.second;
    f_scenario_reduction_config.second = nullptr;
  }
  
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
  
  if (!probsLoaded || poolProbabilities.size() != nbScenarios) {
    // Create uniform weights
    poolProbabilities.assign(nbScenarios, 1.0 / nbScenarios);
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
      // Try to read k (optional - if present, will trigger init_representative_pool)
      ScenarioIndex k = 0;
      bool has_k = false;
      try {
        auto kVar = cfgGroup.getVar("k");
        if (!kVar.isNull()) {
          kVar.getVar(&k);
          has_k = true;
        }
      } catch (...) {
        // k not found, continue without it
      }
      
      // Try to read ell (optional, default 2.0)
      float ell = 2.0f;
      try {
        auto ellVar = cfgGroup.getVar("ell");
        if (!ellVar.isNull()) {
          ellVar.getVar(&ell);
        }
      } catch (...) {
        // ell not found, use default
      }
      
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
      
      // If k was provided, we need to set up for scenario reduction
      if (has_k && k > 0) {
        // Generate missing configs
        if (!block_cfg) {
          block_cfg = generate_default_cfl_config(k, ell);
          if (!block_cfg) {
            throw std::runtime_error("Failed to generate default BlockConfig");
          }
          if (!block_cfg->f_static_variables_Configuration) {
            throw std::runtime_error("generate_default_cfl_config did not set f_static_variables_Configuration");
          }
        }
        if (!solver_cfg) {
          solver_cfg = generate_default_solver_config();
          if (!solver_cfg) {
            throw std::runtime_error("Failed to generate default BlockSolverConfig");
          }
        }
        
        // Set the configuration
        try {
          set_scenario_reduction_config(block_cfg, solver_cfg);
          k_value = k;
          
        } catch (const std::exception& e) {
          // Clean up allocated configs on error
          if (!block_cfg) delete block_cfg;
          if (!solver_cfg) delete solver_cfg;
          throw std::runtime_error(std::string("Failed to set scenario reduction config: ") + e.what());
        }
        
        // Store k for later use (after all data is loaded)
        pending_k_value = k;
      } else if (block_cfg || solver_cfg) {
        // Configs provided without k - just store them for later use
        if (!block_cfg) {
          block_cfg = generate_default_cfl_config(1, ell); // Dummy k=1
        }
        if (!solver_cfg) {
          solver_cfg = generate_default_solver_config();
        }
        set_scenario_reduction_config(block_cfg, solver_cfg);
      }
    }
  } catch (const std::exception& e) {
    // Log warning but continue - scenario reduction config is optional
    std::cerr << "Warning: Failed to load ScenarioReductionConfig: " << e.what() << std::endl;
  }
  
  // After all data is loaded, apply pending scenario reduction if requested
  if (pending_k_value > 0) {
    apply_pending_scenario_reduction();
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
  if (f_scenario_reduction_config.first && f_scenario_reduction_config.second) {
    netCDF::NcGroup cfgGroup = group.addGroup("ScenarioReductionConfig");
    
    // Serialize k if available
    if (k_value > 0) {
      auto kVar = cfgGroup.addVar("k", netCDF::ncInt);
      kVar.putVar(&k_value);
    }
    
    // Serialize ell if available from BlockConfig
    if (f_scenario_reduction_config.first && f_scenario_reduction_config.first->f_static_variables_Configuration) {
      auto* ell_config = dynamic_cast<SimpleConfiguration<double>*>(
          f_scenario_reduction_config.first->f_static_variables_Configuration);
      if (ell_config) {
        float ell = static_cast<float>(ell_config->f_value);
        auto ellVar = cfgGroup.addVar("ell", netCDF::ncFloat);
        ellVar.putVar(&ell);
      }
    }
    
    // Serialize BlockConfig
    auto blockGroup = cfgGroup.addGroup("BlockConfig");
    f_scenario_reduction_config.first->serialize(blockGroup);
    
    // Serialize BlockSolverConfig  
    auto solverGroup = cfgGroup.addGroup("BlockSolverConfig");
    f_scenario_reduction_config.second->serialize(solverGroup);
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
void DiscreteScenarioSet::init_representative_pool( ScenarioIndex k )
{
  // Validate the k parameter
  validate_k_parameter(k);
  
  // Ensure configuration exists (create with defaults if needed)
  ensure_configuration_exists(k);
  
  // Clean up any existing pool
  empty_pool();
  
  // Set up pool parameters
  currentScenarioIndex = 0;
  
  // Apply scenario reduction using the configuration
  apply_scenario_reduction();
  
  // Update pool size based on the selected scenarios
  poolSize = static_cast<ScenarioIndex>(scenarioIndexes.size());
  
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
  if (!f_scenario_reduction_config.first || !f_scenario_reduction_config.second) {
    return false;
  }
  
  try {
    // Try to get the k parameter and validate it
    DiscreteScenarioSet::ScenarioIndex k = get_k_parameter(f_scenario_reduction_config.first);
    
    // If k is valid, scenario reduction can be used
    return (k > 0 && k <= nbScenarios);
  } catch (const std::exception& e) {
    // If parameter extraction fails, don't use scenario reduction
    return false;
  }
}

// Apply scenario reduction
void DiscreteScenarioSet::apply_scenario_reduction()
{
  if (!f_scenario_reduction_config.first || !f_scenario_reduction_config.second) {
    throw std::runtime_error("Missing scenario reduction configuration");
  }
  
  // Extract the k parameter (number of scenarios to select)
  DiscreteScenarioSet::ScenarioIndex k = get_k_parameter(f_scenario_reduction_config.first);
  
  // DEBUG: Log k parameter
  #ifndef NDEBUG
  std::cout << "DEBUG [apply_scenario_reduction]: k (max scenarios to select) = " << k << std::endl;
  std::cout << "DEBUG [apply_scenario_reduction]: total scenarios available = " << nbScenarios << std::endl;
  #endif
  
  // Extract ell parameter for distance calculations from BlockConfig
  float ell = DEFAULT_ELL_VALUE;
  if (f_scenario_reduction_config.first->f_static_variables_Configuration) {
    auto* ell_config = dynamic_cast<SimpleConfiguration<double>*>(
        f_scenario_reduction_config.first->f_static_variables_Configuration);
    if (ell_config) {
      ell = static_cast<float>(ell_config->f_value);
    }
  }
  
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
    auto transport_costs = compute_transport_cost_matrix(n_scenarios, scenario_size, ell);
    
    // Load the CFL problem into the block FIRST
    #ifndef NDEBUG
    std::cout << "DEBUG [apply_scenario_reduction]: Creating CFL block with:" << std::endl;
    std::cout << "  - n_facilities = " << n_scenarios << std::endl;
    std::cout << "  - n_customers = " << n_scenarios << std::endl;
    std::cout << "  - balanced = false" << std::endl;
    std::cout << "  - max_open_facilities = " << k << std::endl;
    #endif
    
    cflBlock->load(
        n_scenarios,       // Number of facilities
        n_scenarios,       // Number of customers
        std::move(capacities),
        std::move(fixed_costs),
        std::move(demands),
        std::move(transport_costs),
        false,             // Not a balanced problem
        k                  // Maximum number of facilities
    );
    
    // Apply the BlockConfig AFTER loading data
    if (f_scenario_reduction_config.first) {
      f_scenario_reduction_config.first->apply(cflBlock.get());
    }
    
    // Generate abstract variables and constraints
    #ifndef NDEBUG
    std::cout << "DEBUG [apply_scenario_reduction]: Generating abstract variables" << std::endl;
    #endif
    cflBlock->generate_abstract_variables();
    
    // Generate constraints with wc=7 to include max facilities constraint
    #ifndef NDEBUG
    std::cout << "DEBUG [apply_scenario_reduction]: Generating abstract constraints with wc=7" << std::endl;
    #endif
    SimpleConfiguration<int> constraint_config(7);  // wc = 7 to generate all constraints
    cflBlock->generate_abstract_constraints(&constraint_config);
    
    // Configure the solver using helper function
    Solver* solver = create_and_configure_solver(cflBlock.get(), ell);
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
      // Extract selected scenarios using helper function
      extract_selected_scenarios(solver, cflBlock.get(), n_scenarios);
      
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


[[nodiscard]] ScenarioGenerator::Scenario DiscreteScenarioSet::get_current_scenario( void )
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

[[nodiscard]] double DiscreteScenarioSet::get_current_scenario_probability( void )
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
DiscreteScenarioSet::get_current_scenario_with_prob()
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
[[nodiscard]] ScenarioGenerator::ScenarioSize DiscreteScenarioSet::get_scenario_size( void )
{
  return scenarioSize;
}

/// Concrete implementation of ScenarioGenerator
DiscreteScenarioSet::DiscreteScenarioSet() { set_seed(DEFAULT_SEED); }

/// Destructor - using RAII principles
DiscreteScenarioSet::~DiscreteScenarioSet() {
  // Reset scenario set
  scenarioSet.resize(boost::extents[0][0]);
  
  // Clear vectors with shrink_to_fit to release memory back to the system
  scenarioIndexes.clear();
  scenarioIndexes.shrink_to_fit();
  
  poolProbabilities.clear();
  poolProbabilities.shrink_to_fit();
  
  // Release the configuration objects
  if (f_scenario_reduction_config.first) {
    delete f_scenario_reduction_config.first;
    f_scenario_reduction_config.first = nullptr;
  }
  
  if (f_scenario_reduction_config.second) {
    delete f_scenario_reduction_config.second;
    f_scenario_reduction_config.second = nullptr;
  }
}

/*--------------------------------------------------------------------------*/
/*-------------------- HELPER METHODS IMPLEMENTATION -----------------------*/
/*--------------------------------------------------------------------------*/


// Validate the k parameter for scenario reduction
void DiscreteScenarioSet::validate_k_parameter(ScenarioIndex k) const
{
  if (k == 0 || k > nbScenarios) {
    throw std::invalid_argument("Invalid k parameter: must be between 1 and " + 
                                std::to_string(nbScenarios));
  }
}

// Ensure scenario reduction configuration exists
void DiscreteScenarioSet::ensure_configuration_exists(ScenarioIndex k)
{
  // If configuration is not initialized, create it with default values
  if (!f_scenario_reduction_config.first || !f_scenario_reduction_config.second) {
    // Create default configurations
    auto* block_cfg = generate_default_cfl_config(k, DEFAULT_ELL_VALUE);
    auto* solver_cfg = generate_default_solver_config("Dupacova");
    
    // Set the configuration
    set_scenario_reduction_config(block_cfg, solver_cfg);
    
    // Verify configuration was created
    if (!f_scenario_reduction_config.first || !f_scenario_reduction_config.second) {
      throw std::runtime_error("Failed to create scenario reduction configuration");
    }
  } else {
    // Configuration exists, but verify it has the correct k parameter
    if (f_scenario_reduction_config.first && f_scenario_reduction_config.first->f_extra_Configuration) {
      auto* k_config = dynamic_cast<SimpleConfiguration<int>*>(f_scenario_reduction_config.first->f_extra_Configuration);
      if (!k_config) {
        // Extra configuration exists but is not SimpleConfiguration<int>, fix it
        delete f_scenario_reduction_config.first->f_extra_Configuration;
        f_scenario_reduction_config.first->f_extra_Configuration = new SimpleConfiguration<int>(k);
      }
    } else if (f_scenario_reduction_config.first) {
      // No extra configuration, add it
      f_scenario_reduction_config.first->f_extra_Configuration = new SimpleConfiguration<int>(k);
    }
  }
  // Case 1: Configuration exists (either just created or from deserialization)
  // TODO: In production, implement a way to update the k parameter in existing
  // BlockConfig. Currently, if configuration was loaded from deserialization,
  // the k parameter from the configuration will be used instead of the one
  // passed to this method. This is a limitation as I don't know if BlockConfig
  // allows this kind of updates.
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
    capacities[i] = 1.0;      // Each facility can serve at most one unit
    fixed_costs[i] = 0.0;     // No fixed cost in scenario reduction
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
  if (f_scenario_reduction_config.second) {
    // Check if this is a MILPSolver configuration
    const auto& solver_names = f_scenario_reduction_config.second->get_SolverNames();
    if (!solver_names.empty()) {
      const std::string& solver_name = solver_names[0];
      bool is_milp_solver = (solver_name == "CPXMILPSolver" || 
                            solver_name == "GRBMILPSolver" || 
                            solver_name == "SCIPMILPSolver" || 
                            solver_name == "HiGHSMILPSolver");
      
      if (is_milp_solver) {
        // For MILPSolver, we need to ensure the ComputeConfig parameters are properly set
        // The configuration should already be properly formatted if loaded from file
        // or created with proper parameter names
      }
    }
    
    f_scenario_reduction_config.second->apply(cflBlock);
  }
  
  // Get the registered solver
  if (cflBlock->get_registered_solvers().empty()) {
    throw std::runtime_error("No solver registered to the block after BlockSolverConfig::apply");
  }
  
  auto* base_solver = cflBlock->get_registered_solvers().front();
  if (!base_solver) {
    throw std::runtime_error("Failed to get solver from block");
  }
  
  // If it's a ScenarioReductionSolver, set the ell parameter
  auto* scenario_solver = dynamic_cast<ScenarioReductionSolver*>(base_solver);
  if (scenario_solver) {
    scenario_solver->set_par(ScenarioReductionSolver::dblEll, static_cast<double>(ell));
  }
  // For MILPSolver, ell parameter is not needed as distances are already in the transport costs
  
  // Return the solver pointer (block owns it)
  return base_solver;
}

// Extract selected scenarios from solver results
void DiscreteScenarioSet::extract_selected_scenarios(const Solver* solver,
                                                    const CapacitatedFacilityLocationBlock* cflBlock,
                                                    ScenarioIndex n_scenarios)
{
  #ifndef NDEBUG
  std::cout << "DEBUG [extract_selected_scenarios]: Starting extraction for solver: " 
            << solver->classname() << std::endl;
  #endif
  
  // Clear existing selection
  scenarioIndexes.clear();
  
  // Check if it's a ScenarioReductionSolver
  auto* scenario_solver = dynamic_cast<const ScenarioReductionSolver*>(solver);
  if (scenario_solver) {
    #ifndef NDEBUG
    std::cout << "DEBUG [extract_selected_scenarios]: Using ScenarioReductionSolver path" << std::endl;
    #endif
    // Get the solution - which scenarios were selected
    const auto& reduced_atoms = scenario_solver->get_reduced_atoms();
    
    // Add the selected scenarios to scenarioIndexes
    for (ScenarioIndex i = 0; i < n_scenarios; ++i) {
      if (reduced_atoms[i]) {
        scenarioIndexes.push_back(i);
      }
    }
  } else {
    #ifndef NDEBUG
    std::cout << "DEBUG [extract_selected_scenarios]: Using MILPSolver path" << std::endl;
    #endif
    
    // For MILPSolver or other solvers, we need to read the y variables from the block
    // First ensure the solver has written the solution to the block
    const_cast<Solver*>(solver)->get_var_solution();
    #ifndef NDEBUG
    std::cout << "DEBUG [extract_selected_scenarios]: Solution written to block" << std::endl;
    #endif
    
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
      std::cout << "DEBUG [extract_selected_scenarios]: y[" << i << "] = " << y_value << std::endl;
      #endif
      
      // Check if facility i is open (y[i] > 0.5)
      if (y_value > 0.5) {
        scenarioIndexes.push_back(i);
      }
    }
  }
  
  #ifndef NDEBUG
  std::cout << "DEBUG [extract_selected_scenarios]: Selected " << scenarioIndexes.size() 
            << " scenarios" << std::endl;
  #endif
  
  // If no scenarios were selected, throw an error
  if (scenarioIndexes.empty()) {
    throw std::runtime_error("No scenarios selected by the reduction algorithm");
  }
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

// Apply pending scenario reduction after deserialization
void DiscreteScenarioSet::apply_pending_scenario_reduction()
{
  if (pending_k_value > 0 && pending_k_value <= nbScenarios) {
    init_representative_pool(pending_k_value);
    pending_k_value = 0;  // Clear the pending value
  }
}


// Generate default BlockConfig for CFL
BlockConfig* DiscreteScenarioSet::generate_default_cfl_config(ScenarioIndex k, float ell) const
{
  auto* config = new BlockConfig(false);  // not differential
  
  // Set k in extra_Configuration
  config->f_extra_Configuration = new SimpleConfiguration<int>(k);
  
  // Set ell in static_variables_Configuration  
  config->f_static_variables_Configuration = new SimpleConfiguration<double>(static_cast<double>(ell));
  
  // Add default constraint generation (all constraints)
  // For now, leave this to the CapacitatedFacilityLocationBlock defaults
  
  
  return config;
}

// Generate default BlockSolverConfig for ScenarioReductionSolver
BlockSolverConfig* DiscreteScenarioSet::generate_default_solver_config(const std::string& algorithm) const
{
  auto* config = new BlockSolverConfig(true);  // differential
  
  // Add ScenarioReductionSolver configuration
  config->add_ComputeConfig("ScenarioReductionSolver", nullptr);
  
  return config;
}

/*--------------------------------------------------------------------------*/
/*------------------ End file DiscreteScenarioSet.cpp ----------------------*/
/*--------------------------------------------------------------------------*/