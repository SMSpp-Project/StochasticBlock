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

#include <Eigen/Dense>
#include "CapacitatedFacilityLocationBlock.h"
#include "ScenarioReductionSolver.h"
#include "Block.h" // For BlockConfig
#include "Solver.h" // For BlockSolverConfig

// Simple configuration class for scenario reduction
/**
 * Simple configuration class for scenario reduction
 * 
 * This avoids API incompatibilities while still providing the necessary functionality.
 */
class ScenarioReductionConfig : public SMSpp_di_unipi_it::Configuration {
public:
    // Simple structure for configuration properties
    struct ConfigProperty {
        std::string name;
        std::variant<int, float, double, std::string> value;
    };
    
    // Nested configuration
    struct NestedConfig {
        std::string type;
        std::vector<ConfigProperty> properties;
    };
    
    // Constructors
    ScenarioReductionConfig(bool diff = true) : Configuration() {}
    
    ScenarioReductionConfig(std::istream& input) : Configuration() {
        load(input);
    }
    
    ScenarioReductionConfig(const ScenarioReductionConfig& old) : Configuration() {
        cfl_config = old.cfl_config;
        solver_config = old.solver_config;
    }
    
    // Required for SMS++ factory
    void load(std::istream& input) override {}
    
    // Create method - required for factory registration
    static SMSpp_di_unipi_it::Configuration* create() { return new ScenarioReductionConfig(); }
    
    // Specialized clone implementation
    ScenarioReductionConfig* clone() const override {
        return new ScenarioReductionConfig(*this);
    }
    
    // Name for the configuration type
    const std::string& private_name() const override {
        static const std::string name = "ScenarioReductionConfig";
        return name;
    }
    
    // Helper to set the k parameter (number of scenarios to select)
    void set_k(int k) {
        cfl_config.properties.push_back({"k", k});
    }
    
    // Helper to set the ell parameter (Wasserstein distance power)
    void set_ell(float ell) {
        cfl_config.properties.push_back({"ell", ell});
    }
    
    // Helper to set the algorithm parameter
    void set_algorithm(const std::string& algorithm) {
        solver_config.properties.push_back({"algorithm", algorithm});
    }
    
    // Custom method to configure based on common parameters
    void configure(int k, float ell = 2.0f, const std::string& algorithm = "Dupacova") {
        // Set CFL configuration
        cfl_config.type = "CFLConfig";
        set_k(k);
        set_ell(ell);
        
        // Set solver configuration
        solver_config.type = "SolverConfig";
        set_algorithm(algorithm);
    }
    
    // Access methods for properties
    int get_k() const {
        for (const auto& prop : cfl_config.properties) {
            if (prop.name == "k" && std::holds_alternative<int>(prop.value)) {
                return std::get<int>(prop.value);
            }
        }
        return 0; // Default
    }
    
    float get_ell() const {
        for (const auto& prop : cfl_config.properties) {
            if (prop.name == "ell" && std::holds_alternative<float>(prop.value)) {
                return std::get<float>(prop.value);
            }
        }
        return 2.0f; // Default
    }
    
    std::string get_algorithm() const {
        for (const auto& prop : solver_config.properties) {
            if (prop.name == "algorithm" && std::holds_alternative<std::string>(prop.value)) {
                return std::get<std::string>(prop.value);
            }
        }
        return "Dupacova"; // Default
    }
    
private:
    // Storage for nested configurations
    NestedConfig cfl_config;
    NestedConfig solver_config;
};

// Registration is handled through the Configuration system itself
// Rather than using the macro directly

// Helper function to create a scenario reduction configuration (static to avoid linker conflicts)
static SMSpp_di_unipi_it::Configuration* create_scenario_reduction_config(int k, float ell = 2.0f, const std::string& algorithm = "Dupacova") {
    // Create and configure ScenarioReductionConfig
    auto config = new ScenarioReductionConfig();
    config->configure(k, ell, algorithm);
    return config;
}

// TestConfig class removed - no longer needed for backward compatibility

// Additional C++ standard library includes not in SMSTypedefs.h
#include <chrono> // C++11: Used for generating unique temporary filenames

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------- FACTORY MANAGEMENT ----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( DiscreteScenarioSet );

/*--------------------------------------------------------------------------*/
/*------------------------- Kmeans clustering ------------------------------*/
/*--------------------------------------------------------------------------*/
/* Naive implementation from scratch of Lloyd's algorithm to solve the
 * k-means optimization problem. There are better-existing libraries that
 * are known to be robust and efficient to solve the k-means optimization
 * problem. The following implementation is simply here as a helper for
 * DiscreteScenarioSet to avoid additional dependencies with external
 * libraries. */

/* Computes the Euclidean distance between two Eigen::VectorXd. */
[[nodiscard]] static double euclideanDistance(const Eigen::VectorXd& vec1,
                                            const Eigen::VectorXd& vec2)
{
  return (vec1 - vec2).norm();
}

/// (Costly) function to find the index of the nearest center to a point
[[nodiscard]] static int nearestCenterIndex(const Eigen::VectorXd& point,
                                           const DiscreteScenarioSet::DiscreteRepresentativePool& centers)
{
  // Handle edge cases
  if (centers.empty()) {
    return -1;  // No centers to compare against
  }
  
  // Use a structured binding to keep track of minimum distance and index
  auto [minDistance, minIndex] = [&centers, &point]() {
    double minDist = std::numeric_limits<double>::max();
    int idx = 0;
    
    for (size_t i = 0; i < centers.size(); i++) {
      // Safety check that vectors are the same size
      if (centers[i].size() != point.size()) {
        continue;  // Skip incompatible vectors
      }
      
      double distance = euclideanDistance(point, centers[i]);
      if (distance < minDist) {
        minDist = distance;
        idx = i;
      }
    }
    
    return std::make_pair(minDist, idx);
  }();
  
  return minIndex;
}

static void kMeans(unsigned int k, DiscreteScenarioSet::PoolMap& pool,
                  DiscreteScenarioSet::DiscreteRepresentativePool& centers,
                  std::vector<int>& labels)
{
  size_t n = pool.size(); // n = nbScenarios
  
  // Safety check for empty pool
  if (n == 0) {
    return;  // Nothing to do with empty pool
  }
  
  unsigned int scenariosize = pool[0].size();
  
  // Ensure centers and labels have the correct size
  if (centers.size() != k) {
    centers.resize(k);
    for (auto& center : centers) {
      if (center.size() != scenariosize) {
        center.resize(scenariosize);
      }
    }
  }
  
  if (labels.size() != n) {
    labels.resize(n, 0);
  }
  
  // Special case: when k=1, just compute the centroid of all points
  if (k == 1) {
    // Set all labels to 0 (there's only one cluster)
    std::fill(labels.begin(), labels.end(), 0);
    
    // Reset the center
    centers[0].setZero();
    
    // Compute the centroid by summing all points
    for (size_t i = 0; i < n; i++) {
      centers[0] += pool[i];
    }
    
    // Divide by the number of points
    centers[0] /= static_cast<double>(n);
    
    // No need for iteration when k=1
    return;
  }
  
  // For k > 1, proceed with standard k-means
  bool changed = true; // Initialize to true to ensure first iteration runs
  int iteration = 0;   // Add iteration counter to avoid infinite loops
  const int MAX_ITERATIONS = 100; // Limit iterations as a safeguard
  
  while (changed && iteration < MAX_ITERATIONS)
  {
    changed = false;
    iteration++;

    // Assign points to the nearest center
    for (size_t i = 0; i < n; i++)
    {
      int newIndex = nearestCenterIndex(pool[i], centers);
      if (newIndex >= 0 && newIndex < static_cast<int>(k) && labels[i] != newIndex)
      {
        labels[i] = newIndex;
        changed = true;
      }
    }

    // Update centers by computing barycenter of each Voronoi cell
    // Use vectors of zeros for each center
    for (auto& center : centers) {
      center.setZero();
    }
    
    // Count points in each cluster
    std::vector<int> counts(k, 0);
    
    // Sum all points in each cluster
    for (size_t i = 0; i < n; i++)
    {
      const int clusterIdx = labels[i];
      if (clusterIdx >= 0 && clusterIdx < static_cast<int>(k)) {
        Eigen::VectorXd& center = centers[clusterIdx];
        
        // Add the point to its cluster center
        center += pool[i];
        
        // Increment count for this cluster
        counts[clusterIdx]++;
      }
    }

    // Compute the average (barycenter) for each cluster
    for (size_t i = 0; i < k; i++)
    {
      // Avoid division by zero
      if (counts[i] > 0) {
        centers[i] /= static_cast<double>(counts[i]);
      }
    }
  }
}

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

void DiscreteScenarioSet::empty_representativePool()
{
  // Only clear the representativePool, not the original probabilities
  // Use a scope-based approach with lambda for better organization
  {
    // We save a copy of the current probabilities
    auto savedProbs = poolProbabilities;
    
    // Clear the pools
    representativePool.clear();
    representativePool.shrink_to_fit();
    poolProbabilities.clear();
    
    // Restore the original probabilities if they exist
    if (!savedProbs.empty()) {
      poolProbabilities = std::move(savedProbs);
    }
  }
  
  // Update the pool type
  if (currentPoolType == PoolType::Continuous) {
    currentPoolType = PoolType::None;
  }
}

[[nodiscard]] bool DiscreteScenarioSet::isempty_representativePool() const
{
  return representativePool.empty();
}

void DiscreteScenarioSet::empty_discretePool()
{
  // Clear the indexes and free memory
  scenarioIndexes.clear();
  scenarioIndexes.shrink_to_fit();
  
  // Update the pool type
  if (currentPoolType == PoolType::Discrete) {
    currentPoolType = PoolType::None;
  }
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
/*--------------------- CLASS DiscreteScenarioSet --------------------------*/
/*--------------------------------------------------------------------------*/

void DiscreteScenarioSet::deserialize( const netCDF::NcGroup & group )
{
  // Reset state to properly handle multiple deserializations
  currentPoolType = PoolType::None;
  currentScenarioIndex = 0;
  poolSize = 0;
  sumPoolWeights = 0.0;

  // Clear existing data
  scenarioSet.resize(boost::extents[0][0]);
  poolProbabilities.clear();
  scenarioIndexes.clear();
  representativePool.clear();
  scenarioReductionConfig.reset(); // Clear any existing configuration
  
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
      try {
        // Try to create a new ScenarioReductionConfig directly
        auto config = new ScenarioReductionConfig();
        
        // Set up configuration properties from the netCDF group
        // Read the 'k' parameter if available
        try {
          netCDF::NcGroup cflGroup = cfgGroup.getGroup("CFLConfig");
          if (!cflGroup.isNull()) {
            netCDF::NcGroupAtt kAtt = cflGroup.getAtt("k");
            if (!kAtt.isNull()) {
              int k_value;
              kAtt.getValues(&k_value);
              config->set_k(k_value);
            }
            
            netCDF::NcGroupAtt ellAtt = cflGroup.getAtt("ell");
            if (!ellAtt.isNull()) {
              float ell_value;
              ellAtt.getValues(&ell_value);
              config->set_ell(ell_value);
            }
          }
        } catch (...) {}
        
        // Read the 'algorithm' parameter if available
        try {
          netCDF::NcGroup solverGroup = cfgGroup.getGroup("SolverConfig");
          if (!solverGroup.isNull()) {
            netCDF::NcGroupAtt algoAtt = solverGroup.getAtt("algorithm");
            if (!algoAtt.isNull()) {
              std::string alg_value;
              algoAtt.getValues(alg_value);
              config->set_algorithm(alg_value);
            }
          }
        } catch (...) {}
        
        // Set the configuration
        scenarioReductionConfig.reset(config);
      } catch (const std::exception& e) {
        // Configuration creation failed, continue without it
        scenarioReductionConfig.reset();
      }
    }
  } catch (const std::exception& e) {
    // Configuration not found or invalid, just continue without it
    scenarioReductionConfig.reset();
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

void DiscreteScenarioSet::init_discrete_pool(ScenarioIndex size)
{
  empty_representativePool();
  sumPoolWeights = 0.0;
  set_poolSize(size);
  currentScenarioIndex = 0;
  currentPoolType = PoolType::Discrete;

  scenarioIndexes.clear();
  
  // We resize instead of reserve to ensure the container has the correct size
  if (size > 0) {
    scenarioIndexes.resize(size);
  }

  // Check if we should use scenario reduction based on configuration
  bool used_reduction = false;
  if (should_use_scenario_reduction(size)) {
    // Try to apply scenario reduction using the configuration
    used_reduction = apply_scenario_reduction(size);
  }
  
  // If scenario reduction was not used or failed, use standard random selection
  if (!used_reduction) {
    generateRandomSubset(nbScenarios, size, scenarioIndexes, rng);
  }

  // Save the total probability weights of the pool in sumPoolWeights
  // Using std::accumulate with lambda for better readability and safety
  sumPoolWeights = 0.0;
  
  if (size > 0) {
    sumPoolWeights = std::accumulate(scenarioIndexes.begin(), scenarioIndexes.end(), 0.0,
      [this](double sum, ScenarioIndex index) -> double {
        return sum + (index < poolProbabilities.size() ? poolProbabilities[index] : 0.0);
      });
  }
}

void DiscreteScenarioSet::init_continuous_pool( ScenarioIndex size )
{
  empty_discretePool();
  set_poolSize( size );
  currentScenarioIndex = 0;
  currentPoolType = PoolType::Continuous;
  
  // Special case: handle size=0
  if (size == 0) {
    // Clear the representative pool and probabilities
    representativePool.clear();
    poolProbabilities.clear();
    return;
  }
  
  poolProbabilities.resize( size , 0 );

  // Viewing every input scenario into an Eigen::VectorXd using a lambda for clarity
  auto createEigenSet = [this]() {
    PoolMap result;
    result.reserve(get_nbScenarios());
    
    for (size_t i = 0; i < get_nbScenarios(); i++) {
      result.emplace_back(Eigen::Map<Eigen::VectorXd>(
        &scenarioSet[i][0], get_scenarioSize()));
    }
    
    return result;
  };
  
  PoolMap eigenSet = createEigenSet();

  // Initialize the representativePool using a random subset of input scenarios
  representativePool.reserve(size);
  std::vector<ScenarioIndex> rand_ind;

  generateRandomSubset(get_nbScenarios(), size, rand_ind, rng);
  
  // Use range-based for loop with structured binding for better readability
  for (const auto& i : rand_ind) {
    representativePool.push_back(Eigen::VectorXd(eigenSet[i]));
  }

  // Lloyd's algo for k-means clustering problem
  // updates in-place representativePool and labels
  std::vector<int> labels(nbScenarios, 0);
  
  // Only run kMeans if we have at least one element
  if (size > 0) {
    // Initialize with uniform probabilities for continuous pool
    std::fill(poolProbabilities.begin(), poolProbabilities.end(), 1.0 / size);
    
    // Only run k-means if we have more than one scenario
    if (nbScenarios > 1) {
      // Save a copy of the original probabilities for later use
      std::vector<double> originalProbs(nbScenarios);
      for (size_t i = 0; i < nbScenarios; ++i) {
        originalProbs[i] = 1.0 / nbScenarios; // Set to uniform
      }
      
      // Run k-means clustering
      kMeans(size, eigenSet, representativePool, labels);
      
      // If k-means successful, compute proper probabilities
      // Reset poolProbabilities first
      std::fill(poolProbabilities.begin(), poolProbabilities.end(), 0.0);
      
      // Count points in each cluster to ensure we have at least one point per cluster
      std::vector<int> cluster_counts(size, 0);
      for (size_t j = 0; j < nbScenarios; j++) {
        if (labels[j] < size) {
          cluster_counts[labels[j]]++;
          poolProbabilities[labels[j]] += originalProbs[j];
        }
      }
      
      // Check if all clusters have at least one point
      bool valid_clustering = true;
      for (size_t i = 0; i < size; i++) {
        if (cluster_counts[i] == 0) {
          valid_clustering = false;
          break;
        }
      }
      
      // If valid clustering, normalize the probabilities
      if (valid_clustering) {
        double sum = 0.0;
        for (size_t i = 0; i < size; i++) {
          sum += poolProbabilities[i];
        }
        
        // If sum is valid, normalize
        if (sum > 0.0) {
          for (size_t i = 0; i < size; i++) {
            poolProbabilities[i] /= sum;
          }
        }
      }
    }
  }
}

[[nodiscard]] ScenarioGenerator::Scenario DiscreteScenarioSet::get_current_scenario( void )
{
  if( currentScenarioIndex >= poolSize )
  {
    throw( std::out_of_range( "Current scenario index is out of range." ) );
  }

  // Use different strategy based on pool type
  switch (currentPoolType) {
    case PoolType::Discrete: {
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
    
    case PoolType::Continuous: {
      // Make sure representativePool has the expected size
      if (representativePool.size() <= currentScenarioIndex)
      {
        throw( std::out_of_range( "Representative pool is too small" ) );
      }
      
      // Transform the currentScenarioIndex-th element of representativePool
      // into a span<const double>
      return Scenario(representativePool[currentScenarioIndex].data(), get_scenario_size());
    }
    
    default:
      throw std::runtime_error("No active pool initialized");
  }
}

[[nodiscard]] double DiscreteScenarioSet::get_current_scenario_probability( void )
{
  if( currentScenarioIndex >= poolSize )
  {
    throw( std::out_of_range( "Current scenario index is out of range." ) );
  }

  // Use different strategy based on pool type
  switch (currentPoolType) {
    case PoolType::Discrete: {
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
    
    case PoolType::Continuous: {
      // Check if currentScenarioIndex is valid for poolProbabilities
      if (currentScenarioIndex >= poolProbabilities.size())
      {
        throw( std::out_of_range( "Probability index is out of range" ) );
      }
    
      return poolProbabilities[currentScenarioIndex];
    }
    
    default:
      throw std::runtime_error("No active pool initialized");
  }
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
    // For discrete pool
    if (currentPoolType == PoolType::Discrete) {
      if (index < scenarioIndexes.size()) {
        const auto scenarioIndex = scenarioIndexes[index];
        if (scenarioIndex < nbScenarios) {
          return Scenario(&scenarioSet[scenarioIndex][0], scenarioSize);
        }
      }
    } 
    // For continuous pool
    else if (currentPoolType == PoolType::Continuous) {
      if (index < representativePool.size()) {
        return Scenario(representativePool[index].data(), scenarioSize);
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
  return currentPoolType != PoolType::None;
}

[[nodiscard]] bool DiscreteScenarioSet::next_scenario( void )
{
  // If poolSize is 0 or no pool is initialized, there are no scenarios to move to
  if (poolSize == 0 || currentPoolType == PoolType::None) {
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

/// Check if scenario reduction should be used based on configuration
bool DiscreteScenarioSet::should_use_scenario_reduction(ScenarioIndex size) const {
  // Skip scenario reduction for zero size
  if (size == 0) return false;
  
  // Skip scenario reduction for size 1, as it's a trivial case
  if (size == 1) return false;
  
  // Check if we have a valid configuration
  if (!scenarioReductionConfig) return false;
  
  // Check if the configuration is a ScenarioReductionConfig
  auto srConfig = dynamic_cast<const ScenarioReductionConfig*>(scenarioReductionConfig.get());
  if (srConfig) {
    // Get the k value from the configuration
    int k = srConfig->get_k();
    
    // If k is specified in the configuration and valid, definitely use scenario reduction
    if (k > 0 && static_cast<ScenarioIndex>(k) <= nbScenarios) {
      return true;
    }
    
    // Even if k is not specified, if we have a valid configuration, it's worth attempting
    // scenario reduction with the provided size parameter
    return true;
  }
  
  // If we have a configuration but it's not a ScenarioReductionConfig, 
  // don't use scenario reduction
  return false;
}

/// Apply scenario reduction based on configuration
bool DiscreteScenarioSet::apply_scenario_reduction(ScenarioIndex size) {
  // Check if we have a valid configuration
  if (!scenarioReductionConfig) {
    return false;
  }

  // Get ell value (power in Wasserstein distance) from configuration
  float ell = 2.0;  // Default to squared Euclidean distance
  
  // Check 'k' parameter (number of scenarios to select) from configuration
  ScenarioIndex k_scenarios = size;
  
  std::string algorithm = "Dupacova"; // Default algorithm
  
  // Try to use the ScenarioReductionConfig
  try {
    auto srConfig = dynamic_cast<ScenarioReductionConfig*>(scenarioReductionConfig.get());
    if (srConfig) {
      // Use the proper methods to access configuration values
      int k_value = srConfig->get_k();
      if (k_value > 0 && static_cast<ScenarioIndex>(k_value) <= nbScenarios) {
        k_scenarios = static_cast<ScenarioIndex>(k_value);
      }
      
      ell = srConfig->get_ell();
      algorithm = srConfig->get_algorithm();
    } else {
      // If the configuration is not of the expected type, use default values
      // This could happen if a different configuration type is provided
      // but we expect specifically ScenarioReductionConfig
    }
  } catch (...) {
    // If anything goes wrong with configuration access, use defaults
  }

  try {    
    #ifdef WITH_CAPACITATED_FACILITY_LOCATION
    // Create a CapacitatedFacilityLocationBlock
    auto cflBlock = std::make_unique<CapacitatedFacilityLocationBlock>();

    // Create the data structures needed for the CFL problem
    const Index n_facilities = nbScenarios;  // Each original scenario is a potential facility
    const Index n_customers = nbScenarios;   // Each original scenario is also a customer

    // Create capacity vector - each facility can serve exactly one customer
    CapacitatedFacilityLocationBlock::DVector capacities(n_facilities, 1.0);
    
    // Create fixed costs vector - all facilities have the same fixed cost
    CapacitatedFacilityLocationBlock::CVector fixed_costs(n_facilities, 0.0);
    
    // Create demands vector - each customer has a demand equal to its probability weight
    CapacitatedFacilityLocationBlock::DVector demands = poolProbabilities;
    
    // Create transportation cost matrix - the cost is the Wasserstein distance between scenarios
    CapacitatedFacilityLocationBlock::CMatrix transport_costs(boost::extents[n_facilities][n_customers]);
    
    // Compute distances between scenarios (transport costs)
    for (ScenarioIndex i = 0; i < n_facilities; ++i) {
      for (ScenarioIndex j = 0; j < n_customers; ++j) {
        // Calculate Euclidean distance between scenarios i and j raised to power ell
        double dist = 0.0;
        for (ScenarioSize d = 0; d < scenarioSize; ++d) {
          double diff = scenarioSet[i][d] - scenarioSet[j][d];
          dist += std::pow(diff, 2);  // Squared difference
        }
        transport_costs[i][j] = std::pow(dist, ell/2.0);  // Convert squared Euclidean to ell-Wasserstein
      }
    }
    
    // Load the CFL problem
    cflBlock->load(n_facilities, n_customers, std::move(capacities), std::move(fixed_costs),
                  std::move(demands), std::move(transport_costs), true, k_scenarios);
    
    // Create and configure the ScenarioReductionSolver
    auto solver = std::make_unique<ScenarioReductionSolver>();
    solver->set_Block(cflBlock.get());
    solver->set_ell(ell);
    
    // Set the algorithm based on what we found in the configuration
    if (algorithm == "Dupacova") {
      solver->set_algorithm(ScenarioReductionSolver::Algorithm::Dupacova);
    } else if (algorithm == "BestFit") {
      solver->set_algorithm(ScenarioReductionSolver::Algorithm::BestFit);
    } else if (algorithm == "FirstFit") {
      solver->set_algorithm(ScenarioReductionSolver::Algorithm::FirstFit);
    } else {
      // Unknown algorithm, fall back to default
      solver->set_algorithm(ScenarioReductionSolver::Algorithm::Dupacova);
    }
    
    // Solve the scenario reduction problem
    int status = solver->compute();
    
    // Check if the solve was successful
    if (status == Solver::kOK) {
      // Get the solution - which scenarios were selected
      const auto& reduced_atoms = solver->get_reduced_atoms();
      
      // Clear the current scenario indexes
      scenarioIndexes.clear();
      
      // Add the selected scenarios to scenarioIndexes
      for (ScenarioIndex i = 0; i < n_facilities; ++i) {
        if (reduced_atoms[i]) {
          scenarioIndexes.push_back(i);
        }
      }
      
      // Calculate sum of weights of selected scenarios for normalization
      sumPoolWeights = 0.0;
      for (const auto& idx : scenarioIndexes) {
        sumPoolWeights += poolProbabilities[idx];
      }
      
      return true;
    } else {
      return false;
    }
    #else // WITH_CAPACITATED_FACILITY_LOCATION not defined
    
    // Fallback to a simpler method: sort scenarios by probability and take the top k
    // This is not scenario reduction in the Wasserstein sense, but provides a deterministic
    // selection that can be useful when the proper solver is not available
    
    // Create a vector of (index, probability) pairs
    std::vector<std::pair<ScenarioIndex, double>> scenario_probs;
    for (ScenarioIndex i = 0; i < nbScenarios; ++i) {
      scenario_probs.emplace_back(i, poolProbabilities[i]);
    }
    
    // Sort by probability in descending order
    std::sort(scenario_probs.begin(), scenario_probs.end(),
              [](const auto& a, const auto& b) { return a.second > b.second; });
    
    // Take the top k_scenarios
    scenarioIndexes.clear();
    sumPoolWeights = 0.0;
    for (ScenarioIndex i = 0; i < std::min(k_scenarios, static_cast<ScenarioIndex>(scenario_probs.size())); ++i) {
      scenarioIndexes.push_back(scenario_probs[i].first);
      sumPoolWeights += scenario_probs[i].second;
    }
    
    return !scenarioIndexes.empty();
    #endif // WITH_CAPACITATED_FACILITY_LOCATION
  } catch (const std::exception& e) {
    // Silently handle errors in scenario reduction and fall back to random selection
    return false;
  }
}

/// Concrete implementation of ScenarioGenerator
DiscreteScenarioSet::DiscreteScenarioSet() { set_seed(DEFAULT_SEED); }

/// Destructor - using RAII principles
DiscreteScenarioSet::~DiscreteScenarioSet() {
  // Clear all containers to free memory
  // Using a lambda to encapsulate the cleanup logic
  auto cleanupContainers = [this]() {
    // Reset scenario set
    scenarioSet.resize(boost::extents[0][0]);
    
    // Clear vectors with shrink_to_fit to release memory back to the system
    scenarioIndexes.clear();
    scenarioIndexes.shrink_to_fit();
    
    representativePool.clear();
    representativePool.shrink_to_fit();
    
    poolProbabilities.clear();
    poolProbabilities.shrink_to_fit();
    
    // Release the configuration
    scenarioReductionConfig.reset();
  };
  
  // Execute cleanup
  cleanupContainers();
}

/*--------------------------------------------------------------------------*/
/*------------------ End file DiscreteScenarioSet.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
