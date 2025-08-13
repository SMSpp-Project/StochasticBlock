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
/*------------------------- STATIC HELPER FUNCTIONS ------------------------*/
/*--------------------------------------------------------------------------*/
/* Static Helper Functions
 * These are internal helper functions not part of the public interface.
 * They provide utility functionality used by the class implementation.
 * 
 */

// Validate the poolSize parameter
static void validate_poolSize_value(DiscreteScenarioSet::ScenarioIndex size, 
                                    DiscreteScenarioSet::ScenarioIndex max_scenarios)
{
  if (size == 0 || size > max_scenarios) {
    throw std::invalid_argument("Invalid pool size parameter: must be between 1 and " + 
                                std::to_string(max_scenarios));
  }
}

// Draw k elements among n with weighted sampling
/* The function generateWeightedRandomSubset draws k elements among n using
 * weighted random sampling without replacement. The chosen indexes are
 * moved into the input variable ind. */
static void generateWeightedRandomSubset( size_t n , size_t k ,
                                          const std::vector< double > & weights,
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

  // Special case: if k == n, return all indices
  if (k == n) {
    ind.resize(n);
    std::iota(ind.begin(), ind.end(), 0);
    return;
  }

  // Use weighted sampling without replacement
  // Create a copy of weights that we can modify
  std::vector<double> working_weights = weights;
  ind.reserve(k);
  
  for (size_t i = 0; i < k; ++i) {
    // Create discrete distribution with current weights
    std::discrete_distribution<ScenarioGenerator::ScenarioIndex> dist(working_weights.begin(), working_weights.end());
    
    // Sample an index
    ScenarioGenerator::ScenarioIndex selected = dist(rng);
    
    // Add to result
    ind.push_back(selected);
    
    // Set weight to 0 to prevent re-selection
    working_weights[selected] = 0.0;
  }
}

// Generate default BlockConfig for CFL (static helper function)
static BlockConfig* generate_default_cfl_config(DiscreteScenarioSet::ScenarioIndex size)
{
  auto* config = new BlockConfig(false);  // not differential
  
  // Set pool size in extra_Configuration
  config->f_extra_Configuration = new SimpleConfiguration<int>(size);
  
  return config;
}

// Create CFL problem data structures
static std::tuple<CapacitatedFacilityLocationBlock::DVector,
                  CapacitatedFacilityLocationBlock::CVector,
                  CapacitatedFacilityLocationBlock::DVector>
create_cfl_problem_data(DiscreteScenarioSet::ScenarioIndex n_scenarios,
                       const std::vector<double>& poolWeights)
{
  CapacitatedFacilityLocationBlock::DVector capacities(n_scenarios);
  CapacitatedFacilityLocationBlock::CVector fixed_costs(n_scenarios);
  CapacitatedFacilityLocationBlock::DVector demands(n_scenarios);
  
  // Set up the basic CFL parameters
  for (DiscreteScenarioSet::ScenarioIndex i = 0; i < n_scenarios; ++i) {
    capacities[i] = 1.0;      // Capacities = weights of the reduced distribution
    fixed_costs[i] = 0.0;     // No fixed cost in scenario reduction
    demands[i] = poolWeights[i];  // Demand equals weight
  }
  
  return std::make_tuple(std::move(capacities), std::move(fixed_costs), std::move(demands));
}

// Extract selected scenarios from CFL block results
static void extract_scenarios_from_cfl_block(const CapacitatedFacilityLocationBlock* cflBlock,
                                             DiscreteScenarioSet::ScenarioIndex n_scenarios,
                                             std::vector<DiscreteScenarioSet::ScenarioIndex>& scenarioIndexes)
{
  #ifndef NDEBUG
  std::cout << "DEBUG [extract_scenarios_from_cfl_block]: Starting extraction from CFL block" << std::endl;
  #endif
  
  // Clear existing selection
  scenarioIndexes.clear();
    
  // Read variable values directly from the block variables
  for (DiscreteScenarioSet::ScenarioIndex i = 0; i < n_scenarios; ++i) {
    // Get the y variable for facility i
    const auto* y_var = cflBlock->get_y(i);
    if (!y_var) {
      throw std::runtime_error("Failed to get y variable for facility " + std::to_string(i));
    }
    
    // Get value from the variable directly
    double y_value = y_var->get_value();
    #ifndef NDEBUG
    std::cout << "DEBUG [extract_scenarios_from_cfl_block]: y[" << i << "] = " << y_value << std::endl;
    #endif
    
    // Check if facility i is open (y[i] > 0.5)
    if (y_value > 0.5) {
      scenarioIndexes.push_back(i);
    }
  }
  
  #ifndef NDEBUG
  std::cout << "DEBUG [extract_scenarios_from_cfl_block]: Selected " << scenarioIndexes.size() 
            << " scenarios: ";
  for (auto idx : scenarioIndexes) std::cout << idx << " ";
  std::cout << std::endl;
  #endif
  
  // If no scenarios were selected, throw an error
  if (scenarioIndexes.empty()) {
    throw std::runtime_error("No scenarios selected by the reduction algorithm");
  }
}

/*--------------------------------------------------------------------------*/
/*----------------------- GETTERS AND SETTERS ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

// Get the indices of currently selected scenarios
const std::vector<DiscreteScenarioSet::ScenarioIndex>& 
DiscreteScenarioSet::get_selected_scenarios() const
{
  if (!is_initialized) {
    throw std::runtime_error("Pool has not been initialized. Call init_random_pool() or init_representative_pool() first.");
  }
  return scenarioIndexes;
}


/*--------------------------------------------------------------------------*/
/*-------------------- HELPER METHODS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

void DiscreteScenarioSet::empty_pool()
{
  // Clear the indexes and free memory
  scenarioIndexes.clear();
  scenarioIndexes.shrink_to_fit();
  
  // Clear normalized weights and free memory
  normalizedPoolWeights.clear();
  normalizedPoolWeights.shrink_to_fit();
  
  // Reset pool state (poolSize is configuration, not actual size)
  currentScenarioIndex = 0;
  // poolSize is NOT reset - it's the configured/desired size
  sumPoolWeights = 0.0;
  is_initialized = false;
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
void DiscreteScenarioSet::set_config(BlockConfig* block_config, BlockSolverConfig* solver_config, ScenarioIndex pool_size)
{
  // Validate and set pool size
  validate_poolSize_value(pool_size, nbScenarios);
  poolSize = pool_size;
  
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
        // Extract poolSize from the BlockConfig's extra configuration
        if (!block_config->f_extra_Configuration) {
          throw std::runtime_error("No extra configuration found containing poolSize parameter");
        }
        
        auto* simple_config = dynamic_cast<SimpleConfiguration<int>*>(block_config->f_extra_Configuration);
        if (!simple_config) {
          throw std::runtime_error("Extra configuration is not a SimpleConfiguration<int>");
        }
        
        ScenarioIndex pool_size = simple_config->f_value;
        
        // Validate the extracted poolSize
        validate_poolSize_value(pool_size, nbScenarios);
        
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
  poolWeights.clear();
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
  bool probsLoaded = ::deserialize(group, "poolWeights", nbScenarios, poolWeights);
  
  if (!probsLoaded) {
    // No probabilities in file, create uniform weights
    poolWeights.assign(nbScenarios, 1.0 / nbScenarios);
  } else if (poolWeights.size() != nbScenarios) {
    // Probabilities were loaded but have wrong size - this is an error
    throw std::invalid_argument("poolWeights size (" + std::to_string(poolWeights.size()) + 
                                ") does not match NumberScenarios (" + std::to_string(nbScenarios) + ")");
  }
  
  // Validate probabilities sum to approximately 1.0
  double sum = std::accumulate(poolWeights.begin(), poolWeights.end(), 0.0);
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
        // No BlockSolverConfig
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
  auto probsVar = group.addVar("poolWeights", netCDF::NcDouble(), nbScenariosDim);
  probsVar.putVar(poolWeights.data());
  
  // Serialize scenario reduction configuration if it exists
  // We serialize if we have at least poolSize or BlockConfig or BlockSolverConfig
  if (poolSize > 0 || f_block_config || f_solver_config) {
    netCDF::NcGroup cfgGroup = group.addGroup("ScenarioReductionConfig");
    
    // Serialize poolSize if available
    if (poolSize > 0) {
      auto poolSizeVar = cfgGroup.addVar("poolSize", netCDF::ncInt);
      poolSizeVar.putVar(&poolSize);
    }
    
    // Serialize ell from internal variable
    auto ellVar = cfgGroup.addVar("ell", netCDF::ncFloat);
    ellVar.putVar(&ell);
    
    // Serialize BlockConfig if it exists
    if (f_block_config) {
      auto blockGroup = cfgGroup.addGroup("BlockConfig");
      f_block_config->serialize(blockGroup);
    }
    
    // Serialize BlockSolverConfig if it exists
    if (f_solver_config) {
      auto solverGroup = cfgGroup.addGroup("BlockSolverConfig");
      f_solver_config->serialize(solverGroup);
    }
  }
}


// Initialize a pool with randomly selected scenarios
void DiscreteScenarioSet::init_random_pool(ScenarioIndex pool_size)
{
  // Clean up any existing pool
  empty_pool();
  
  // Validate and set pool size
  validate_poolSize_value(pool_size, nbScenarios);
  poolSize = pool_size;
  currentScenarioIndex = 0;
  
  // Generate random indices for the pool using weighted sampling
  if (pool_size > 0) {
    generateWeightedRandomSubset(nbScenarios, pool_size, poolWeights, scenarioIndexes, rng);
    
    update_pool_weights();
      
    // Mark the pool as initialized
    is_initialized = true;
  }
}

// Initialize a representative pool using scenario reduction
void DiscreteScenarioSet::init_representative_pool( ScenarioIndex target_pool_size )
{
  // Step 1: Validate the target_pool_size parameter
  validate_poolSize_value(target_pool_size, nbScenarios);
  
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
    auto [capacities, fixed_costs, demands] = create_cfl_problem_data(n_scenarios, poolWeights);
    
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
        target_pool_size   // Maximum number of facilities to open
    );
    
    // Step 3.1: Apply BlockConfig if it exists, else generate minimal default
    if (f_block_config) {
      // Trust that the user provided appropriate BlockConfig for CFL
      f_block_config->apply(cflBlock.get());
    } else {
      // Generate minimal default config (mostly empty, just target_pool_size in extra_Configuration)
      auto* default_config = generate_default_cfl_config(target_pool_size);
      default_config->apply(cflBlock.get());
      delete default_config;
    }
    
    // Step 3.2: Generate abstract variables (to read the solution)
    cflBlock->generate_abstract_variables();
    
    // Apply BlockSolverConfig to register the solver with the Block
    // IMPORTANT: apply() transfers ownership of the BlockSolverConfig to the solver
    f_solver_config->apply(cflBlock.get());
    f_solver_config.release();  // Release our ownership - the solver now owns it
    
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
    
    // Ensure solver solution is written to CFL Block variables
    solver->get_var_solution();
    
    // Extract the solution from the CFL block
    extract_scenarios_from_cfl_block(cflBlock.get(), n_scenarios, scenarioIndexes);
  } else {
    // No solver config - use baseline method (select top target_pool_size by weight)
    // BlockConfig alone is not useful for baseline selection
    apply_baseline_selection(target_pool_size);
  }
  
  // Update pool weights and finalize
  update_pool_weights();
  
  // Reset current index (poolSize remains as configured)
  currentScenarioIndex = 0;
  
  // Mark the pool as initialized if we have scenarios
  if (!scenarioIndexes.empty()) {
    is_initialized = true;
  }
}



[[nodiscard]] bool DiscreteScenarioSet::next_scenario( void )
{
  // If no pool is initialized or empty, there are no scenarios to move to
  if (!is_initialized || scenarioIndexes.empty()) {
    return false;
  }
  
  if (currentScenarioIndex < scenarioIndexes.size() - 1)
  {
    // Use prefix increment for efficiency
    ++currentScenarioIndex;
    return true; // Successfully moved to the next scenario
  }
  return false; // No more scenario in scenarioPool to move to
}


/// Concrete implementation of ScenarioGenerator
DiscreteScenarioSet::DiscreteScenarioSet() { 
  set_seed(1337);  // Default seed for reproducibility
}

/// Destructor
DiscreteScenarioSet::~DiscreteScenarioSet() {
  // Reset scenario set
  scenarioSet.resize(boost::extents[0][0]);
  
  // Clear vectors with shrink_to_fit to release memory back to the system
  scenarioIndexes.clear();
  scenarioIndexes.shrink_to_fit();
  
  poolWeights.clear();
  poolWeights.shrink_to_fit();
  
  normalizedPoolWeights.clear();
  normalizedPoolWeights.shrink_to_fit();
  
  // Release the configuration objects following SMS++ ownership patterns
  
  // BlockConfig: We own our clone, safe to delete
  if (f_block_config) {
    delete f_block_config;
    f_block_config = nullptr;
  }
  
  // BlockSolverConfig: unique_ptr handles cleanup
  // If ownership was transferred via release(), this will be nullptr
}

/*--------------------------------------------------------------------------*/
/*-------------------- HELPER METHODS IMPLEMENTATION -----------------------*/
/*--------------------------------------------------------------------------*/


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


// Apply baseline selection method
void DiscreteScenarioSet::apply_baseline_selection(ScenarioIndex target_size)
{
  // Create pairs of (index, weight) for sorting
  std::vector<std::pair<ScenarioIndex, double>> indexed_weights;
  indexed_weights.reserve(nbScenarios);
  
  for (ScenarioIndex i = 0; i < nbScenarios; ++i) {
    double weight = (i < poolWeights.size()) ? poolWeights[i] : 1.0 / nbScenarios;
    indexed_weights.emplace_back(i, weight);
  }
  
  // Sort by weight in descending order (highest weights first)
  std::sort(indexed_weights.begin(), indexed_weights.end(),
            [](const auto& a, const auto& b) { return a.second > b.second; });
  
  // Select the top target_size scenarios
  scenarioIndexes.clear();
  scenarioIndexes.reserve(target_size);
  
  for (ScenarioIndex i = 0; i < target_size && i < nbScenarios; ++i) {
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
    sumPoolWeights += poolWeights[idx];
  }
  
  // Populate normalized weights
  normalizedPoolWeights.clear();
  normalizedPoolWeights.reserve(scenarioIndexes.size());
  
  if (sumPoolWeights > 0.0) {
    for (const auto& idx : scenarioIndexes) {
      normalizedPoolWeights.push_back(poolWeights[idx] / sumPoolWeights);
    }
  } else {
    // Fallback to uniform if all weights are zero
    double uniform_weight = 1.0 / scenarioIndexes.size();
    normalizedPoolWeights.assign(scenarioIndexes.size(), uniform_weight);
  }
}

/*--------------------------------------------------------------------------*/
/*------------------ End file DiscreteScenarioSet.cpp ----------------------*/
/*--------------------------------------------------------------------------*/