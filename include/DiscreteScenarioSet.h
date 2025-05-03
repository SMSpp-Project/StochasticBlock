/*--------------------------------------------------------------------------*/
/*-------------------- File DiscreteScenarioSet.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class DiscreteScenarioSet that is an
 * implementation of ScenarioGenerator suited to the case where the input
 * distribution is contained in a netCDF file as a collection of vectors.
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
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __DiscreteScenarioSet
 #define __DiscreteScenarioSet
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ScenarioGenerator.h"

#include <Eigen/Dense>
#include <random>
#include <span>
#include <optional>
#include <variant>
#include <algorithm>
#include <numeric>
#include <utility>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/// User-defined literal for probability percentages (must be at namespace or global scope)
/** This literal allows writing probabilities as percentages, e.g., 25.0_pct */
constexpr double operator"" _pct(long double percentage) {
    return static_cast<double>(percentage / 100.0);
}
/*--------------------------------------------------------------------------*/
/*--------------------- CLASS DiscreteScenarioSet --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/

/// DiscreteScenarioSet to sample from a collection of scenarios
/** DiscreteScenarioSet is an implementation of the ScenarioGenerator class.
 * As such, it gives methods to sample from an input distribution and manipulate
 * a scenarioPool.
 *
 * In the specific context of DiscreteScenarioSet, the distribution to sample
 * from is assumed to be a discrete probability distribution characterized
 * by a collection of scenarios. Scenarios are assumed to be contained in a
 * netCDF file, and DiscreteScenarioSet gives method to deserialize the scenarios
 * from the netCDF file. The deserialized scenarios are stored in a
 * boost::multi_array< double, 2 >.
 *
 * DiscreteScenarioSet considers that the (deserialized) pool can be handled
 * by two distinct approaches for scenario reduction: the first draws a *subset* 
 * of the input scenarioPool (discrete approach), while the second *constructs* 
 * a set of representative scenarios from the input (continuous approach). Several 
 * functions thus have two different behaviors, whether the pool is initialized by 
 * the first or second method.
 *
 * The method get_current_scenario() allows the user to query one element in
 * the pool, whether the first or second method was used.
 *
 * The first method (discrete approach) outputs scenarios that are part of the 
 * input, so the pool is simply characterized by the set of indexes of the drawn 
 * scenarios. See init_discrete_pool(...).
 *
 * The second method (continuous approach) outputs scenarios that can be different 
 * from the ones that were inputted. By default, DiscreteScenarioSet uses a native
 * implementation of k-means clustering. In that case, another container is
 * created to hold the constructed scenarios. See init_continuous_pool(...). */

class DiscreteScenarioSet : public ScenarioGenerator
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:
/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *  @{ */

 /// Container for the deserialized scenario pool
 /** Every scenario is assumed to have the same dimension. As the number of
  * scenarios becomes known whenever we deserialize the data, the scenario
  * pool is of known size at this point. Hence, the choice to store it inside
  * a boost::multi_array. See the function init_random_pool()*/
 using DiscreteScenarioPool = boost::multi_array< double , 2 >;

 /// Container for the representative pool
 /** The representative pool is made of Eigen::VectorXd for ease of linear
  * algebra manipulations. See the function init_representative_pool(). */
 using DiscreteRepresentativePool = std::vector< Eigen::VectorXd >;

 /// type Point for Eigen::VectorXd
 /** For ease of linear algebra manipulations, we use Eigen::VectorXd to
  * hold a scenario. As scenarios are assumed to be present in memory
  * already (as boost::multi_array< double, 2 >), they are converted into
  * Eigen::VectorXd using Eigen::Map. */
 using Point = Eigen::Map< Eigen::VectorXd >;

 /// Container for the centers and the points to be clustered
 /** Lightweight container for a collection of scenarios represented as
  * Eigen::VectorXd. */
 using PoolMap = std::vector< Point >;
 
 /// Type to represent a scenario with its probability
 using ScenarioWithProbability = std::pair<Scenario, double>;
 

/** @} ---------------------------------------------------------------------*/
/*----------- CONSTRUCTING AND DESTRUCTING DiscreteScenarioSet -------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing DiscreteScenarioSet
 *  @{ */

 DiscreteScenarioSet();

 /// deserialize a discrete distribution from a netCDF group
 /** Implementation of the "third-level" pure virtual function deserialize of
  * ScenarioGenerator.h. Assumes that there is a two-dimensional variable
  *  \p Scenario contained inside a netCDF NcGroup. One dimension NumberScenarios
  * corresponds to the number of input scenarios characterizing the input discrete
  * probability distribution. The second dimension ScenarioSize is the dimension
  * of the Euclidean space (R^d) representing a single scenario. We
  * deserialize the scenarios into a boost::multi_array< double, 2 > as the two
  * dimensions become known once the file has been read.
  *
  * The scenarios are associated with another variable ScenarioProbabilities.
  * If ScenarioProbabilities is present in the group, then it is saved in a
  * std::vector< double > called scenarioProbabilities. If ScenarioProbabilities
  * is *not* present in the group, then uniform weights are assumed, that is,
  * scenarioProbabilities is a vector of size nbScenarios where each component
  * is equal to 1.0 / nbScenarios. */
 void deserialize( const netCDF::NcGroup & group ) override;

 virtual ~DiscreteScenarioSet();

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS INHERITED FROM ScenarioGenerator.h ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Functions inherited from ScenarioGenerator.h
 *  @{ */

 void set_seed( unsigned long seed ) override;

 /// Function to select a discrete subset from the input scenario pool
 /** The function init_discrete_pool selects a subset of scenarios among the
  * ones that were deserialized from the input. It saves this subset of indices
  * in the variable scenarioIndexes.
  * 
  * This approach preserves the original scenarios without generating new ones,
  * making it appropriate for scenario reduction when you need to maintain the
  * original scenarios and just want to select a representative subset.
  */
 void init_discrete_pool( ScenarioIndex sampleSize ) override;

 /// Function to compute a continuous approximation for scenario representation
 /** The function init_continuous_pool( ScenarioIndex size ) computes a
  * set of scenarios that approximates the input set of scenarios according
  * to some optimization criterion or statistical properties.
  *
  * By default, in DiscreteScenarioSet a native (and naive) implementation of
  * k-means clustering is used. While k-means aims to minimize variance within 
  * clusters (not necessarily Wasserstein distance), it serves as a general 
  * method for creating representative scenarios.
  *
  * Scenario reduction techniques, which specifically aim to minimize Wasserstein
  * distance between original and reduced distributions, can also be implemented
  * through this interface.
  *
  * Kmeans splits the input scenarios into size clusters. Each cluster has a 
  * (bary)center, and each scenario is associated with a label, which is its 
  * nearest center.
  *
  * Knowing the cluster centers and the labels, the representativePool is
  * made of the centers. The std::vector< double > poolProbabilities is such
  * that each component p_i is equal to
  *  p_i = sum( input_weights_with_label_equal_to_i ).
  * 
  * This method typically creates new scenarios that weren't in the original set
  * but are constructed to be good statistical representatives.
  */
 void init_continuous_pool( ScenarioIndex sampleSize ) override;

 /// Function for retrieving the current scenario.
 /** Checks that the internal variable currentScenarioIndex is within bounds,
  * then converts the currentScenarioIndex-th row of the scenario pool as a
  * Scenario, that is as a std::span< const double >.
  *
  * If a discrete pool has been initialized (via init_discrete_pool), we output 
  * the scenario as a span of the scenarioIndex[currentScenarioIndex]-th row of the
  * scenarioSet.
  *
  * If a continuous pool has been initialized (via init_continuous_pool), we output 
  * the scenario as a span of the currentScenarioIndex-th component of the 
  * representativePool. */
 Scenario get_current_scenario( void ) override;

 /// Function to query the probability weight of the current scenario
 /** When sampling a pool, what are the weights of the drawn scenarios?
  * When using get_scenario_probabilities, we return "the" probability weight
  * inside the pool and not the input probability weight.
  * This applies for both discrete and continuous pool approaches.
  *
  * Choices made:
  * 1) For discrete pools (init_discrete_pool): Take the input weight and 
  * normalize it, that is given a scenario with an input_weight, we compute 
  * its new_weight by:
  *  "new_weight = input_weight / sum( input_weights_in_the_pool )".
  * 
  * 2) For continuous pools (init_continuous_pool): For k-means clustering,
  * the "pool weight" is the proportion of points affected by the center.
  *
  * In the discrete case, after scaling, we return the currentScenarioIndex-th
  * input_weight from the deserialized data in scenarioProbabilities.
  * In the continuous case, we return the currentScenarioIndex-th element of the
  * std::vector< double > poolProbabilities that is constructed with the
  * representative scenarios.
  * */
 double get_current_scenario_probability( void ) override;

 /// Move currentScenarioIndex to the next scenario
 /** Whether the pool has been created via init_discrete_pool(...)
  * or via init_continuous_pool(...), the function next_scenario() 
  * behaves the same: it increments by 1 the currentScenarioIndex if 
  * there is still a scenario left in the pool and return true, 
  * otherwise it returns false. */
 bool next_scenario( void ) override;

 /// return the dimension of the scenarios
 /** Every scenario (vector in some Euclidean space R^d) is assumed to have
  * the same dimension. The dimension has been saved in scenarioSize when
  * deserializing the input discrete distribution. */
 ScenarioSize get_scenario_size( void ) override;

/** @} ---------------------------------------------------------------------*/
/*--------------------- GETTERS  FOR PRIVATE FIELDS ------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Getters for some private fields
 *  @{ */

 /// get a reference to nbScenarios
 const ScenarioIndex & get_nbScenarios() const;

 /// get a reference to scenarioSize
 const ScenarioSize & get_scenarioSize() const;
 
 /// Get current scenario with its probability as a pair
 [[nodiscard]] ScenarioWithProbability get_current_scenario_with_prob();
 
 /// Try to get a scenario by index, returns nullopt if index is invalid
 [[nodiscard]] std::optional<Scenario> try_get_scenario(ScenarioIndex index) const;
 
 /// Check if the scenario pool has been initialized
 [[nodiscard]] bool is_pool_initialized() const;

/** @} ---------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Private fields
 *  @{ */

 /// Current index in the pool
 ScenarioIndex currentScenarioIndex{0};

 /// Number of different scenarios in the scenario pool
 ScenarioIndex nbScenarios;

 /// Size of a scenario
 ScenarioSize scenarioSize;

 /// Container for Scenario-s
 DiscreteScenarioPool scenarioSet;

 /// Pool size
 /** The variable poolSize is initialized when
  * init_representative_pool( size_t size ) or init_rando_pool( size_t size ) is
  * used. */
 ScenarioIndex poolSize = 0;

 /// Random generator
 std::mt19937 rng;
 
 /// Enum for pool type
 enum class PoolType { 
     None,      ///< No pool initialized
     Discrete,  ///< Discrete pool
     Continuous ///< Continuous pool
 };
 
 /// Current pool type
 PoolType currentPoolType{PoolType::None};
 
 /// Compile-time constants
 static constexpr double DEFAULT_EPSILON = 1e-10;
 static constexpr unsigned long DEFAULT_SEED = 1337;

/** @} ---------------------------------------------------------------------*/
/*--------------------- FIELDS FOR REPRESENTATIVE POOL ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Fields for the representative pool
   * When the scenario pool is made of representative scenarios that
   * are in general different from the input scenarios, the scenario pool
   * has its own container of type DiscreteRepresentativePool with
   * associated vector of probability weights poolProbabilities.
 * @{ */

 /// Container for representative pool of scenarios
 DiscreteRepresentativePool representativePool;

 /// Probabilities of scenarios inside the pool
 /** No matter the chosen method to make the pool, poolProbabilities is
  * the vector of associated probability weights. It is different from
  * scenarioProbabilities which is the vector of probability weights
  * of the input scenarios, which is in general a different set than the
  * pool. */
 std::vector< double > poolProbabilities;

/** @} ---------------------------------------------------------------------*/
/*------------------------- FIELDS FOR DISCRETE POOL -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Fields for the discrete pool
   * When the scenario pool is made of a discrete subset of the input scenarios,
   * it is characterized by a std::vector< ScenarioIndex > and the probability
   * weights can be deduced from the input weights saved in
   * scenarioProbabilities and the sumPoolWeights of the scenarios that
   * belong to the scenario pool.
 * @{ */

 /// holder for the sum of the weights inside the discrete pool
 /** Variable which holds the sum of the weights of the scenarios that were
  * chosen to be part of the discrete pool, see init_discrete_pool(...).
  * This variable is set back to 0.0 if the continuous pool is used,
  * see init_continuous_pool(...). */
 double sumPoolWeights;

 /// Indexes of the discrete pool
 std::vector< ScenarioIndex > scenarioIndexes;

/** @} ---------------------------------------------------------------------*/
/*-------------------- HELPER METHODS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name helper methods of the class
 * Miscellaneous functions
 * @{ */

 /// update the variable poolSize
 /** Whenever a size for the pool has been given, update the
  * variable poolSize accordingly. Also ensures that the desired
  * poolSize is possible, that is, it should be less or equal than the
  * number of input scenarios.
  *
  * We normalize the weights of each scenario in the pool so that their
  * weights sum up to one. That is, we have
  *  new_weight = input_weight / sum( input_weights_of_pool_scenarios ),
  * where input_weight refers to the weight of a given scenario contained in
  * the vector scenarioProbabilities.
  *
  * So for the use of get_current_scenario_probability(), we
  * save in memory sumPoolWeight, equal to the sum of the input weights of
  * the scenarios chosen to be part of the pool. */
 void set_poolSize( ScenarioIndex size );

 /// empty the representativePool of scenarios
 /** Function to clear the internal representativePool. */
 void empty_representativePool();

 /// "empty" the discrete pool
 /** As the discrete pool of DiscreteScenarioSet is simply a subset of indices
  * scenarioIndexes, we clear scenarioIndexes. The variable sumPoolWeights
  * is initialized to 0.0 */
 void empty_discretePool();

 /// function to check if the representativePool of scenarios is empty
 bool isempty_representativePool() const;

  SMSpp_insert_in_factory_h;

};   // end( class DiscreteScenarioSet )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 } // end( namespace SMSpp_di_unipi_it )

#endif /* __DiscreteScenarioSet */

/*--------------------------------------------------------------------------*/
/*------------------- End file DiscreteScenarioSet.h -----------------------*/
/*--------------------------------------------------------------------------*/

