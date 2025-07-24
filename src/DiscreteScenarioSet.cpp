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
#include "ScenarioReductionSolver.h"
#include "Block.h" 
#include "Solver.h"
#include "BlockSolverConfig.h" // For BlockSolverConfig class
#include "Configuration.h" // For Configuration class

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
  // Clean up existing configurations
  if (f_scenario_reduction_config.first) {
    delete f_scenario_reduction_config.first;
  }
  
  if (f_scenario_reduction_config.second) {
    delete f_scenario_reduction_config.second;
  }
  
  // Set the new configurations
  f_scenario_reduction_config.first = block_config;
  f_scenario_reduction_config.second = solver_config;
}

// Set scenario reduction configuration with k parameter
void DiscreteScenarioSet::set_scenario_reduction_config(BlockConfig* block_config, BlockSolverConfig* solver_config, ScenarioIndex k)
{
  // First set k_value with validation
  set_k_value(k);
  
  // Then set the configurations
  set_scenario_reduction_config(block_config, solver_config);
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
      // Deserialize BlockConfig and BlockSolverConfig from their respective subgroups
      BlockConfig* block_cfg = nullptr;
      BlockSolverConfig* solver_cfg = nullptr;
      
      try {
        auto blockGroup = cfgGroup.getGroup("BlockConfig");
        if (!blockGroup.isNull()) {
          block_cfg = dynamic_cast<BlockConfig*>(Configuration::new_Configuration(blockGroup));
        }
      } catch (...) {
        // BlockConfig group not found, continue
      }
      
      try {
        auto solverGroup = cfgGroup.getGroup("BlockSolverConfig");
        if (!solverGroup.isNull()) {
          solver_cfg = dynamic_cast<BlockSolverConfig*>(Configuration::new_Configuration(solverGroup));
        }
      } catch (...) {
        // BlockSolverConfig group not found, continue
      }
      
      // If we got valid configurations, store them
      if (block_cfg && solver_cfg) {
        set_scenario_reduction_config(block_cfg, solver_cfg);
        
        // Sync k_value with the configuration
        try {
          k_value = get_k_parameter(block_cfg);
        } catch (const std::exception& e) {
          // If k parameter extraction fails, keep default k_value = 0
          k_value = 0;
        }
      } else {
        // Clean up partial configurations
        delete block_cfg;
        delete solver_cfg;
      }
    }
  } catch (const std::exception& e) {
    // Configuration group not found or invalid, continue without it
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
    
    // Apply the BlockConfig to the CapacitatedFacilityLocationBlock
    if (f_scenario_reduction_config.first) {
      f_scenario_reduction_config.first->apply(cflBlock.get());
    }
    
    // Set up the scenario selection problem parameters
    ScenarioIndex n_scenarios = nbScenarios;
    ScenarioSize scenario_size = scenarioSize;
    
    // Create CFL problem data using helper function
    auto [capacities, fixed_costs, demands] = create_cfl_problem_data(n_scenarios);
    
    // Compute the transport cost matrix using helper function
    auto transport_costs = compute_transport_cost_matrix(n_scenarios, scenario_size, ell);
    
    // Load the CFL problem into the block
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
    
    // Configure the solver using helper function
    ScenarioReductionSolver* solver = create_and_configure_solver(cflBlock.get(), ell);
    
    // Solve the scenario reduction problem
    int status = solver->compute();
    
    // Check if the solve was successful
    if (status == Solver::kOK) {
      // Extract selected scenarios using helper function
      extract_selected_scenarios(solver, n_scenarios);
      
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

void DiscreteScenarioSet::create_scenario_reduction_config(
    ScenarioIndex k,
    float ell,
    const std::string& algorithm,
    double rho,
    bool shuffle,
    unsigned long random_seed)
{
  // Check if configuration is already set
  if (f_scenario_reduction_config.first != nullptr || 
      f_scenario_reduction_config.second != nullptr) {
    // Configuration already exists, do nothing
    return;
  }
  
  // Validate mandatory k parameter
  if (k == 0) {
    throw std::invalid_argument("k must be positive");
  }
  
  // Validate optional parameters
  if (ell <= 0.0f) {
    throw std::invalid_argument("ell must be positive");
  }
  
  if (algorithm.empty()) {
    throw std::invalid_argument("algorithm cannot be empty");
  }
  
  // Note: k_value member variable is set during deserialization or via set_k_value()
  // The configuration created here stores k in f_extra_Configuration for consistency
  
  // Create temporary txt files for the configurations
  std::string block_config_content = R"(# BlockConfig for scenario reduction
BlockConfig

0  # not a differential configuration

# static constraints Configuration
* # [none]

# dynamic constraints Configuration
* # [none]

# static variables Configuration - ell parameter
SimpleConfiguration<double>
)" + std::to_string(ell) + R"(

# dynamic variables Configuration
* # [none]

# objective Configuration
* # [none]

# is_feasible Configuration
* # [none]

# is_optimal Configuration
* # [none]

# solution Configuration
* # [none]

# extra Configuration - k parameter for scenario reduction
SimpleConfiguration<int>
)" + std::to_string(k) + R"(
)";

  std::string solver_config_content = R"(# BlockSolverConfig for scenario reduction
BlockSolverConfig

0  # not a differential configuration

1  # number of Solver
ScenarioReductionSolver

1  # number of ComputeConfig

# ComputeConfig for ScenarioReductionSolver
ComputeConfig

1  # differential mode

0  # number of integer parameters

1  # number of double parameters
dblRho )" + std::to_string(rho) + R"(

1  # number of string parameters
strAlgorithm )" + algorithm + R"(

0  # number of vector-of-int parameters

0  # number of vector-of-double parameters

0  # number of vector-of-string parameters

# extra Configuration with shuffle and random_seed
SimpleConfiguration<std::pair<int,int>>
)" + std::to_string(shuffle ? 1 : 0) + " " + std::to_string(random_seed) + R"(
)";

  // Create temporary files
  std::string block_config_file = "/tmp/sr_block_config_" + 
      std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".txt";
  std::string solver_config_file = "/tmp/sr_solver_config_" + 
      std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".txt";
  
  // Write configurations to files
  std::ofstream block_file(block_config_file);
  block_file << block_config_content;
  block_file.close();
  
  std::ofstream solver_file(solver_config_file);
  solver_file << solver_config_content;
  solver_file.close();
  
  // Load configurations from files
  BlockConfig* block_config = dynamic_cast<BlockConfig*>(
      Configuration::deserialize(block_config_file));
  BlockSolverConfig* solver_config = dynamic_cast<BlockSolverConfig*>(
      Configuration::deserialize(solver_config_file));
  
  // Clean up temporary files
  std::remove(block_config_file.c_str());
  std::remove(solver_config_file.c_str());
  
  // Store the configuration pair
  f_scenario_reduction_config = std::make_pair(block_config, solver_config);
}

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
  // Case 2: If configuration is not initialized, create it with default values
  if (!f_scenario_reduction_config.first || !f_scenario_reduction_config.second) {
    // Create configuration with k and default values for other parameters
    create_scenario_reduction_config(k);
    
    // Verify configuration was created
    if (!f_scenario_reduction_config.first || !f_scenario_reduction_config.second) {
      throw std::runtime_error("Failed to create scenario reduction configuration");
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
  Eigen::VectorXd diff = scenario1 - scenario2;
  double norm;
  
  // Euclidean norm for ell=2 (more efficient)
  if (ell == 2.0f) {
    norm = diff.norm();
  } else {
    // For other values of ell, compute the ell-norm
    norm = std::pow((diff.array().abs().pow(ell)).sum(), 1.0/ell);
  }
  
  // Transportation cost = ell-power of the norm
  return std::pow(norm, ell);
}

// Create and configure the scenario reduction solver
ScenarioReductionSolver*
DiscreteScenarioSet::create_and_configure_solver(CapacitatedFacilityLocationBlock* cflBlock,
                                                 float ell) const
{
  // Apply the BlockSolverConfig to register and configure the solver
  if (f_scenario_reduction_config.second) {
    f_scenario_reduction_config.second->apply(cflBlock);
  }
  
  // Get the registered solver (ScenarioReductionSolver)
  auto* base_solver = cflBlock->get_registered_solvers().front();
  auto* solver = dynamic_cast<ScenarioReductionSolver*>(base_solver);
  
  if (!solver) {
    throw std::runtime_error("Failed to get ScenarioReductionSolver from block");
  }
  
  // Set the ell parameter (needed for distance calculations)
  solver->set_par(ScenarioReductionSolver::dblEll, static_cast<double>(ell));
  
  // Return the solver pointer (block owns it)
  return solver;
}

// Extract selected scenarios from solver results
void DiscreteScenarioSet::extract_selected_scenarios(const ScenarioReductionSolver* solver,
                                                    ScenarioIndex n_scenarios)
{
  // Get the solution - which scenarios were selected
  const auto& reduced_atoms = solver->get_reduced_atoms();
  
  // Clear existing selection
  scenarioIndexes.clear();
  
  // Add the selected scenarios to scenarioIndexes
  for (ScenarioIndex i = 0; i < n_scenarios; ++i) {
    if (reduced_atoms[i]) {
      scenarioIndexes.push_back(i);
    }
  }
  
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

/*--------------------------------------------------------------------------*/
/*------------------ End file DiscreteScenarioSet.cpp ----------------------*/
/*--------------------------------------------------------------------------*/