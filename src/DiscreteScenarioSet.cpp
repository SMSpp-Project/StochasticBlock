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
  // Use a structured binding to keep track of minimum distance and index
  auto [minDistance, minIndex] = [&centers, &point]() {
    double minDist = std::numeric_limits<double>::max();
    int idx = 0;
    
    for (size_t i = 0; i < centers.size(); i++) {
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
  unsigned int scenariosize = pool[0].size();
  
  bool changed;
  do
  {
    changed = false;

    // Assign points to the nearest center
    for (int i = 0; i < n; i++)
    {
      int newIndex = nearestCenterIndex(pool[i], centers);
      if (labels[i] != newIndex)
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
    for (int i = 0; i < n; i++)
    {
      const int clusterIdx = labels[i];
      Eigen::VectorXd& center = centers[clusterIdx];
      
      // Add the point to its cluster center
      center += pool[i];
      
      // Increment count for this cluster
      counts[clusterIdx]++;
    }

    // Compute the average (barycenter) for each cluster
    for (int i = 0; i < k; i++)
    {
      // Avoid division by zero
      if (counts[i] > 0) {
        centers[i] /= static_cast<double>(counts[i]);
      }
    }
  } while (changed);
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

  generateRandomSubset(nbScenarios, size, scenarioIndexes, rng);

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
  };
  
  // Execute cleanup
  cleanupContainers();
}

/*--------------------------------------------------------------------------*/
/*------------------ End file DiscreteScenarioSet.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
