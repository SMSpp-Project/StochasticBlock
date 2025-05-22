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

#include "DiscreteScenarioSet.h"  // Already includes Eigen/Dense

// Include required components for scenario reduction
#include "CapacitatedFacilityLocationBlock.h"
#include "ScenarioReductionSolver.h"
#include "Block.h"  // For BlockConfig
#include "Solver.h" // For solver class
#include "BlockSolverConfig.h" // For BlockSolverConfig class
#include "MILPSolver.h" // For MILPSolver class
#include "Configuration.h" // For Configuration class

using namespace SMSpp_di_unipi_it;

// Needed for unique temporary filename generation and file operations
#include <chrono>
#include <cstdio>  // for std::remove

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

// Get block configuration for scenario reduction
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
  
  // Validate the k parameter
  if (k_value == 0) {
    throw std::runtime_error("k parameter must be set before using scenario reduction. "
                              "Use set_k_value() or set_scenario_reduction_config() with k parameter.");
  }
  
  if (k_value > nbScenarios) {
    throw std::runtime_error("Invalid k parameter: cannot exceed number of scenarios");
  }
  
  return k_value;
}

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS DiscreteScenarioSet --------------------------*/
/*--------------------------------------------------------------------------*/

void DiscreteScenarioSet::deserialize( const netCDF::NcGroup & group )
{
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
  
  // Compute the two dimensions of the scenarioPool
  ::deserialize_dim( group , "NumberScenarios" , nbScenarios , false );
  ::deserialize_dim( group , "ScenarioSize" , scenarioSize , false );

  // Deserialize the Scenarios inside the scenarioPool
  scenarioSet.resize( boost::extents[ nbScenarios ][ scenarioSize ] );
  ::deserialize( group , "Scenarios" , scenarioSet , true , false );

  // If weights are not present, assume uniform weights
  // Use a lambda to create uniform weights when needed
  auto createUniformWeights = [this]() {
    // Create vector with uniform weights
    return std::vector<double>(nbScenarios, 1.0 / nbScenarios);
  };

  // Load probabilities or create uniform ones
  bool probsLoaded = ::deserialize(group, "poolProbabilities", nbScenarios, poolProbabilities);
  
  if (!probsLoaded || poolProbabilities.size() != nbScenarios) {
    poolProbabilities = createUniformWeights();
  }
  
  // Check for scenario reduction configuration group
  try {
    // Check if the ScenarioReductionConfig group exists
    netCDF::NcGroup cfgGroup = group.getGroup("ScenarioReductionConfig");
    if (!cfgGroup.isNull()) {
      // Create configuration objects
      BlockConfig* block_cfg = new BlockConfig();
      BlockSolverConfig* solver_cfg = new BlockSolverConfig();
      
      // Extract parameters directly from attributes in the ScenarioReductionConfig group
      
      // Extract k parameter
      try {
        netCDF::NcGroupAtt kAtt = cfgGroup.getAtt("k");
        if (!kAtt.isNull()) {
          kAtt.getValues(&k_value);
        }
      } catch (...) {}
      
      // These parameters are stored internally in the netCDF attributes
      // In a production environment, the values would come from these attributes
      // For our implementation, we'll use the BlockConfig and BlockSolverConfig
      // objects as containers, but we won't be using set_parameter/get_parameter
      
      // Store the configurations
      set_scenario_reduction_config(block_cfg, solver_cfg);
    }
  } catch (const std::exception& e) {
    // Configuration not found or invalid, just continue without it
  }
}

void DiscreteScenarioSet::serialize(const netCDF::NcGroup& group) const
{
  // First serialize the base data
  // (Scenarios and poolProbabilities are assumed to be serialized elsewhere)
  
  // Serialize the scenario reduction configuration if it exists
  if (f_scenario_reduction_config.first && f_scenario_reduction_config.second) {
    netCDF::NcGroup cfgGroup = group.addGroup("ScenarioReductionConfig");
        
    // k parameter - number of scenarios to select
    cfgGroup.putAtt("k", netCDF::NcType::nc_INT, k_value);
    
    // TODO: This serialize method needs to be reworked to properly extract
    // the other parameters (ell, algorithm, rho, shuffle, random_seed) from
    // the BlockConfig and BlockSolverConfig objects and serialize them.
    // Currently these parameters are stored in the configuration objects
    // but not serialized to the netCDF file.
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
void DiscreteScenarioSet::init_random_pool(size_t pool_size)
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
void DiscreteScenarioSet::init_representative_pool()
{
  // Check if we have a valid configuration
  if (!f_scenario_reduction_config.first || !f_scenario_reduction_config.second) {
    throw std::runtime_error("No scenario reduction configuration available");
  }
  
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
  
  // TODO: Extract these parameters from the BlockConfig and BlockSolverConfig
  // For now, use default values
  float ell = DEFAULT_ELL_VALUE;
  std::string algo = "Dupacova";
  double rho = DEFAULT_RHO_VALUE;
  bool shuffle = false;
  unsigned long random_seed = DEFAULT_SEED;
  
  // In a real implementation, we would extract from the configuration objects
  // Try to extract from BlockConfig if available
  if (f_scenario_reduction_config.first) {
    // TODO: Extract ell parameter from BlockConfig's static variables Configuration
  }
  
  // Try to extract from BlockSolverConfig if available
  if (f_scenario_reduction_config.second) {
    // TODO: Extract algorithm, rho, shuffle, random_seed from BlockSolverConfig's ComputeConfig
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
    
    // Create CFL problem parameters
    CapacitatedFacilityLocationBlock::DVector capacities(n_scenarios);
    CapacitatedFacilityLocationBlock::CVector fixed_costs(n_scenarios);
    CapacitatedFacilityLocationBlock::DVector demands(n_scenarios);
    CapacitatedFacilityLocationBlock::CMatrix transport_costs(boost::extents[n_scenarios][n_scenarios]);
    
    // Map scenarios to Eigen matrix for easier distance calculations
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> 
        all_scenarios(scenarioSet.data(), n_scenarios, scenario_size);
    
    // Set up the basic CFL parameters
    for (ScenarioIndex i = 0; i < n_scenarios; ++i) {
        capacities[i] = 1.0;      // Each facility can serve at most one unit
        fixed_costs[i] = 0.0;     // No fixed cost in scenario reduction
        demands[i] = poolProbabilities[i];  // Demand equals probability weight
    }
    
    // Compute the distance/transportation cost matrix
    for (ScenarioIndex i = 0; i < n_scenarios; ++i) {
        for (ScenarioIndex j = 0; j < n_scenarios; ++j) {
            if (i == j) {
                transport_costs[i][j] = 0.0;
                continue;
            }
            
            // Calculate p-norm distance between scenarios
            Eigen::VectorXd diff = all_scenarios.row(i) - all_scenarios.row(j);
            double norm;
            
            // Euclidean norm for ell=2 (more efficient)
            if (ell == 2.0f) {
                norm = diff.norm();
            } else {
                // For other values of ell, compute the ell-norm
                norm = std::pow((diff.array().abs().pow(2)).sum(), 1.0/2.0);
            }
            
            // Transportation cost = ell-power of the norm
            transport_costs[i][j] = std::pow(norm, ell);
        }
    }
    
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
    
    // Create and configure the scenario reduction solver
    auto solver = std::make_unique<ScenarioReductionSolver>();
    solver->set_Block(cflBlock.get());
    solver->set_ell(ell);
    
    // Configure the solver based on the selected algorithm
    if (algo == "MILP") {
      // Create a MILP solver for exact scenario selection
      auto milpSolver = std::make_unique<MILPSolver>();
      
      // Register the MILP solver with the block
      milpSolver->set_Block(cflBlock.get());
      
      // Solve using the MILP solver
      if (milpSolver->compute() != Solver::kOK) {
        throw std::runtime_error("MILP solver computation failed");
      }
      
      // Extract the selected scenarios
      CapacitatedFacilityLocationBlock::CntSolution y(n_scenarios);
      cflBlock->get_facility_solution(y.begin());
      
      // Add the selected scenarios to scenarioIndexes
      for (ScenarioIndex i = 0; i < n_scenarios; ++i) {
        if (y[i] > 0.5) { // Using threshold to account for numerical precision
          scenarioIndexes.push_back(i);
        }
      }
    } else {
      // Configure the heuristic solver based on the algorithm
      if (algo == "Dupacova") {
        solver->set_algorithm(ScenarioReductionSolver::Algorithm::Dupacova);
      } else if (algo == "BestFit") {
        solver->set_algorithm(ScenarioReductionSolver::Algorithm::BestFit);
      } else if (algo == "FirstFit") {
        solver->set_algorithm(ScenarioReductionSolver::Algorithm::FirstFit);
      } else {
        // Unknown algorithm, fall back to default
        solver->set_algorithm(ScenarioReductionSolver::Algorithm::Dupacova);
      }
      
      // Apply additional configuration parameters
      solver->set_rho(rho);
      solver->set_shuffle(shuffle);
      solver->set_random_seed(static_cast<unsigned int>(random_seed));
      
      // Solve the scenario reduction problem
      int status = solver->compute();
      
      // Check if the solve was successful
      if (status == Solver::kOK) {
        // Get the solution - which scenarios were selected
        const auto& reduced_atoms = solver->get_reduced_atoms();
        
        // Add the selected scenarios to scenarioIndexes
        for (ScenarioIndex i = 0; i < n_scenarios; ++i) {
          if (reduced_atoms[i]) {
            scenarioIndexes.push_back(i);
          }
        }
      } else {
        throw std::runtime_error("Scenario reduction solver failed with status: " + std::to_string(status));
      }
    }
    
    // If no scenarios were selected, throw an error
    if (scenarioIndexes.empty()) {
      throw std::runtime_error("No scenarios selected by the reduction algorithm");
    }
    
    // Calculate sum of weights of selected scenarios for normalization
    sumPoolWeights = 0.0;
    for (const auto& idx : scenarioIndexes) {
      sumPoolWeights += poolProbabilities[idx];
    }
    
  } catch (const std::exception& e) {
    // If scenario reduction fails, throw a runtime error
    throw std::runtime_error("Scenario reduction failed: " + std::string(e.what()));
  }
}

void DiscreteScenarioSet::init_discrete_pool(ScenarioIndex sampleSize)
{
  // Since we've removed the continuous pool approach, init_discrete_pool 
  // now has the same functionality as init_random_pool but with additional
  // support for scenario reduction when available

  // Clean up existing pool
  empty_pool();
  
  // Initialize basic pool parameters
  sumPoolWeights = 0.0;
  set_poolSize(sampleSize);
  currentScenarioIndex = 0;

  // Check if we should use scenario reduction based on configuration
  bool used_reduction = false;
  if (should_use_scenario_reduction(sampleSize)) {
    try {
      // Try to apply scenario reduction using the configuration
      apply_scenario_reduction();
      used_reduction = true;
    } catch (const std::exception& e) {
      // If scenario reduction fails, fall back to random selection
      std::cerr << "Warning: Scenario reduction failed: " << e.what() << std::endl;
      std::cerr << "Falling back to random scenario selection." << std::endl;
      used_reduction = false;
    }
  }
  
  // If scenario reduction was not used or failed, use standard random selection
  if (!used_reduction) {
    // We resize instead of reserve to ensure the container has the correct size
    if (sampleSize > 0) {
      scenarioIndexes.resize(sampleSize);
    }
    generateRandomSubset(nbScenarios, sampleSize, scenarioIndexes, rng);
  }

  // Save the total probability weights of the pool in sumPoolWeights
  // Using std::accumulate with lambda for better readability and safety
  sumPoolWeights = 0.0;
  
  if (sampleSize > 0 && !scenarioIndexes.empty()) {
    sumPoolWeights = std::accumulate(scenarioIndexes.begin(), scenarioIndexes.end(), 0.0,
      [this](double sum, ScenarioIndex index) -> double {
        return sum + (index < poolProbabilities.size() ? poolProbabilities[index] : 0.0);
      });
    
    // Mark the pool as initialized
    is_initialized = true;
  }
}

// This method is removed as it's obsolete - we no longer support continuous pools
void DiscreteScenarioSet::init_continuous_pool( ScenarioIndex sampleSize )
{
  // For backward compatibility, this method now redirects to init_random_pool
  // to maintain the same interface but avoid the continuous pool approach
  init_random_pool(sampleSize);
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
  
  // Update k_value in the class (the only field that remains as a member)
  k_value = k;
  
  // Create temporary txt files for the configurations
  std::string block_config_content = R"(# BlockConfig for scenario reduction
BlockConfig

0  # not a differential configuration

# static constraints Configuration - k parameter for scenario reduction
SimpleConfiguration<int>
)" + std::to_string(k) + R"(

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

# extra Configuration
* # [none]
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

/*--------------------------------------------------------------------------*/
/*------------------ End file DiscreteScenarioSet.cpp ----------------------*/
/*--------------------------------------------------------------------------*/