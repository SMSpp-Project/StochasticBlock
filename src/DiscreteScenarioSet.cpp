/*--------------------------------------------------------------------------*/
/*-------------------- File DiscreteScenarioSet.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Class DiscreteScenarioSet that is an implementation of ScenarioGenerator
 * suited to the case where the input distribution is contained in a netCDF file
 * as a collection of vectors. Associated with the header file
 * DiscreteScenarioSet.h
 *
 * \author Antonio Frangioni \n Dipartimento di Informatica \n Universita' di
 *         Pisa \n
 *
 * \author Benoit Tran \n Dipartimento di Informatica \n Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
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
/*------------------------- Kmeans clustering ------------------------------*/
/*--------------------------------------------------------------------------*/
/* Naive implementation from scratch of Lloyd's algorithm to solve the
 * kmeans optimization problem. There are better existing librairies that
 * are known to be robust and efficient to solve the kmeans optimization
 * problem. The following implementation is simply here as a helper for
 * DiscreteScenarioSet to avoid additional dependencies with external
 * librairires. */

/* Computes the euclidean distance between two Eigen::VectorXd. */
static double euclideanDistance(const Eigen::VectorXd &vec1,
                                const Eigen::VectorXd &vec2)
{
  return (vec1 - vec2).norm();
}

/// (Costly) function to find the index of the nearest center to a point
static int nearestCenterIndex(const Eigen::VectorXd &point,
                      DiscreteScenarioSet::DiscreteRepresentativePool &centers)
{
  double minDistance = std::numeric_limits<double>::max();
  int index = 0;
  for (size_t i = 0; i < centers.size(); i++)
  {
    double distance = euclideanDistance(point, centers[i]);
    if (distance < minDistance)
    {
      minDistance = distance;
      index = i;
    }
  }
  return index;
}

static void kMeans(unsigned int k, DiscreteScenarioSet::PoolMap &pool,
                   DiscreteScenarioSet::DiscreteRepresentativePool &centers,
                   std::vector<int> &labels)
{
  size_t n = pool.size(); // n = nbScenarios
  unsigned int scenariosize = pool[0].size();
  std::vector<double> counts;
  counts.resize(k);

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
    for (auto &center : centers)
    {
      center.setZero();
    }

    for (int i = 0; i < n; i++)
    {
      Eigen::VectorXd &center = centers[labels[i]];

      for (size_t j = 0; j < scenariosize; j++)
      {
        center[j] += pool[i][j];
      }
      counts[labels[i]]++;
    }

    for (int i = 0; i < k; i++)
    {
      for (size_t j = 0; j < scenariosize; j++)
      {
        centers[i][j] /= counts[i];
      }
    }
  } while (changed);
}

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS DiscreteScenarioSet --------------------------*/
/*--------------------------------------------------------------------------*/

void DiscreteScenarioSet::deserialize(const netCDF::NcGroup &group)
{
  // Compute the two dimensions of the scenarioPool
  ::deserialize_dim(group, "NumberScenarios",
                    nbScenarios, false);

  ::deserialize_dim(group, "ScenarioSize",
                    scenarioSize, false);

  // Deserialize the Scenarios inside the scenarioPool
  scenarioSet.resize(boost::extents[nbScenarios][scenarioSize]);
  ::deserialize(group, "Scenarios", scenarioSet, true,
                false);

  // If weights are not present, assume uniform weights
  if (!::deserialize(group, "ScenarioProbabilities",
                     nbScenarios, scenarioProbabilities))
  {
    scenarioProbabilities.resize(nbScenarios, 1.0 / nbScenarios);
  }
  // maybe instead consider empty if not given and change accordingly
}

// Implementation for setting the seed of the pseudo-random number generator
void DiscreteScenarioSet::set_seed(unsigned long seed)
{
  rng.seed(seed);
}

// Draw k elements among n
/* The function generateRandomSubset draws k elements among n by use of
 * the std::shuffle function and the internal rng. The chosen indexes are
 * moved into the input variable ind. */
static void generateRandomSubset(size_t n, size_t k,
                            std::vector<ScenarioGenerator::ScenarioIndex> &ind,
                            std::mt19937 &rng)
{
  if (k > n)
  {
    throw std::invalid_argument("k must be less or equal than n");
  }
  // elem = [1,2, ..., n]
  std::vector<ScenarioGenerator::ScenarioIndex> elem(n);
  std::iota(elem.begin(), elem.end(), 0);

  // Shuffle elem randomly using our rng
  std::shuffle(elem.begin(), elem.end(), rng);

  // Move the first k ind into indexes
  std::move(elem.begin(), elem.begin() + k, std::back_inserter(ind));
}

void DiscreteScenarioSet::update_poolSize(ScenarioIndex size)
{
  if (size > nbScenarios)
  {
    throw std::out_of_range("The desired sample size is greater than "
                            "the number of available number of different "
                            "scenarios");
  }
  poolSize = size;
}

void DiscreteScenarioSet::empty_representativePool()
{
  representativePool.clear();
  representativePool.shrink_to_fit();
}

void DiscreteScenarioSet::empty_randomPool()
{
  scenarioIndexes.clear();
  scenarioIndexes.shrink_to_fit();
  sumPoolWeights = 0.0;
}

void DiscreteScenarioSet::init_random_pool(ScenarioIndex size)
{
  empty_representativePool();
  update_poolSize(size);

  // the pool will be made of a subset of the input scenarios
  scenarioIndexes.resize(size);

  // Update in-place scenarioIndexes
  generateRandomSubset(nbScenarios, size, scenarioIndexes, rng);

  // Save the total probability weights of the pool in sumPoolWeights
  sumPoolWeights = 0.0;
  for (auto i{0}; i < size; i++)
  {
    sumPoolWeights += scenarioProbabilities[scenarioIndexes[i]];
  }

  currentScenarioIndex = 0;
}

void DiscreteScenarioSet::init_representative_pool(ScenarioIndex size)
{
  empty_randomPool();
  update_poolSize(size);
  poolProbabilities.resize(size);

  // Convert every input scenario into an Eigen::VectorXd
  PoolMap eigenSet;
  eigenSet.reserve(nbScenarios);
  for (size_t i = 0; i < scenarioSet.shape()[0]; i++)
  {
    eigenSet.emplace_back(Eigen::Map<Eigen::VectorXd>(&scenarioSet[i][0],
                                                      scenarioSize));
  }

  // Initialize the representativePool, using a random subset of input scenari
  representativePool.reserve(size);
  std::vector<ScenarioIndex> rand_ind;

  generateRandomSubset(nbScenarios, size, rand_ind, rng);
  for (auto i : rand_ind)
  {
    representativePool.push_back(Eigen::VectorXd(eigenSet[i]));
  }

  // Lloyd's algo for kmeans clustering problem
  // updates in-place representativePool and labels
  std::vector<int> labels(nbScenarios, 0);
  kMeans(size, eigenSet, representativePool, labels);

  // compute the poolProbabilities from the labels and input weights
  for (size_t j{0}; j < nbScenarios; j++)
  {
    poolProbabilities[labels[j]] += scenarioProbabilities[j];
  }

  currentScenarioIndex = 0;
}

/// Bool function to check is the representative set is empty
bool DiscreteScenarioSet::isempty_rep_set()
{
  return representativePool.empty();
}

ScenarioGenerator::Scenario DiscreteScenarioSet::get_current_scenario(void)
{
  if (currentScenarioIndex >= poolSize)
  {
    throw std::out_of_range("Current scenario index is out of range.");
  }

  if (isempty_rep_set())
  {
    // transform the scenarioIndexes[currentScenarioIndex]-th row of
    // scenarioSet into a span< const double >
    return Scenario(&scenarioSet[scenarioIndexes[currentScenarioIndex]][0],
                    get_scenario_size());
  }
  else
  {
    // transform the currentScenarioIndex-th element of representativePool
    // into a span< const double >
    return Scenario(representativePool[currentScenarioIndex].data(),
                    get_scenario_size());
  }
}

double DiscreteScenarioSet::get_current_scenario_probability(void)
{
  if (currentScenarioIndex >= poolSize)
  {
    throw std::out_of_range("Current scenario index is out of range.");
  }
  if (isempty_rep_set())
  {
    return scenarioProbabilities[currentScenarioIndex] / sumPoolWeights;
  }
  else
  {
    return poolProbabilities[currentScenarioIndex];
  }
}

bool DiscreteScenarioSet::next_scenario(void)
{
  if (currentScenarioIndex < poolSize - 1)
  {
    currentScenarioIndex++;
    return true; // Successfully moved to the next scenario
  }
  return false; // No more scenario in scenarioPool to move to
}

/// Implementation for retrieving the size of a scenario
ScenarioGenerator::ScenarioSize DiscreteScenarioSet::get_scenario_size(void)
{
  return scenarioSize;
}

/// Concrete implementation of ScenarioGenerator
DiscreteScenarioSet::DiscreteScenarioSet()
{
  set_seed(1337);
}

/// Destructor
DiscreteScenarioSet::~DiscreteScenarioSet()
{
}

/** @} ---------------------------------------------------------------------*/
/*-------------------------- FACTORY MANAGEMENT ----------------------------*/
/*--------------------------------------------------------------------------*/
SMSpp_insert_in_factory_cpp_0(DiscreteScenarioSet);

/*--------------------------------------------------------------------------*/
/*------------------ End file DiscreteScenarioSet.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
