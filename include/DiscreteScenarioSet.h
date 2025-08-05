/*--------------------------------------------------------------------------*/
/*-------------------- File DiscreteScenarioSet.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class DiscreteScenarioSet that is an
 * implementation of ScenarioGenerator suited to the case where the input
 * distribution is contained in a netCDF file as a collection of vectors.
 * 
 * The class provides methods for scenario selection and management:
 * - Scenario Pool Selection:
 *   - Using random selection (always available)
 *   - Using scenario reduction via Wasserstein distance minimization: 
 *     - Forward selection method (Dupačová)
 *     - Local search algorithms: FirstFit and BestFit
 *     - MILP-based optimization approach
 *
 * The scenario reduction functionality can be configured through a Configuration 
 * object loaded from a netCDF file or set programmatically. It integrates with 
 * the CapacitatedFacilityLocationBlock module to implement scenario selection methods 
 * based on the Wasserstein distance metric.
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
#include "CapacitatedFacilityLocationBlock.h"
#include "ScenarioReductionSolver.h"
#include "Configuration.h"
#include "Block.h"           // For BlockConfig
#include "BlockSolverConfig.h" // For BlockSolverConfig
#include "Solver.h"          // For Solver base class

#include <Eigen/Dense>  // For vector operations

// Additional C++ standard library includes not in SMSTypedefs.h
#include <random>    // For random number generation (std::mt19937)
#include <optional>  // For representing optional values
#include <variant>   // For type-safe unions

// <span> is included in ScenarioGenerator.h
// <numeric> and <utility> are included in SMSTypedefs.h

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
 * As such, it provides methods to sample from an input distribution and manipulate
 * a scenarioPool.
 *
 * In the specific context of DiscreteScenarioSet, the distribution to sample
 * from is assumed to be a discrete probability distribution characterized
 * by a collection of scenarios. Scenarios are assumed to be contained in a
 * netCDF file, and DiscreteScenarioSet provides methods to deserialize the scenarios
 * from the netCDF file. The deserialized scenarios are stored in a
 * boost::multi_array< double, 2 >.
 *
 * DiscreteScenarioSet implements scenario selection from the input pool:
 * 
 * - Selects a *subset* of the input scenarioPool via:
 *   - Random selection (always available) via init_random_pool()
 *   - Scenario reduction via init_representative_pool() (requires proper configuration)
 *   - The subset is characterized by the set of indexes of the selected scenarios
 *
 * The method get_current_scenario() allows the user to query one element in
 * the selected pool. The get_current_scenario_probability() method returns
 * the normalized probability of the current scenario.
 */

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
  * a boost::multi_array. */
 using DiscreteScenarioPool = boost::multi_array< double , 2 >;

 /// Type to hold a scenario reference (used within implementation)
 /** For ease of linear algebra manipulations, we use Eigen::VectorXd when
  * scenario operations need to be performed. */
 using Point = Eigen::Map< Eigen::VectorXd >;
 
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
  * is equal to 1.0 / nbScenarios. 
  * 
  * If a "ScenarioReductionConfig" group is found during deserialization,
  * it will be processed as follows:
  * - k (int): Optional. If provided and > 0, init_representative_pool(k) 
  *   will be called after deserialization completes.
  * - ell (float): Optional. Power for Wasserstein distance (default 2.0).
  * - BlockConfig/: Optional. Configuration for CapacitatedFacilityLocationBlock.
  *   If k is provided but BlockConfig is not, a default will be generated.
  * - BlockSolverConfig/: Optional. Configuration for the solver.
  *   If k is provided but BlockSolverConfig is not, a default ScenarioReductionSolver 
  *   configuration will be generated.
  * 
  * Note: If k is provided, scenario reduction will be automatically applied
  * after the scenarios are loaded.
  */
 void deserialize( const netCDF::NcGroup & group ) override;

 /**
  * @brief Serialize the object to a netCDF group
  * 
  * Extends the base class serialization to include scenario reduction configuration.
  * Creates a "ScenarioReductionConfig" group if the configuration exists.
  * 
  * This method:
  * 1. First calls the base class serialize to save basic scenario data
  * 2. Checks if scenario reduction configuration exists
  * 3. If it exists, creates a "ScenarioReductionConfig" group
  * 4. Creates "BlockConfig" and "SolverConfig" subgroups
  * 5. Writes parameters like k, ell, and algorithm to these groups
  * 
  * Implementation notes:
  * - Only serializes the configuration if it exists
  * - Follows SMS++ netCDF serialization patterns
  * - Preserves all configuration parameters
  * 
  * @param group The netCDF group to serialize the object to
  */
 void serialize(netCDF::NcGroup& group) const override;

 virtual ~DiscreteScenarioSet();

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS INHERITED FROM ScenarioGenerator.h ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Functions inherited from ScenarioGenerator.h
 *  @{ */

 void set_seed( unsigned long seed ) override;

 /// Function for retrieving the current scenario.
 /** Checks that the internal variable currentScenarioIndex is within bounds,
  * then converts the currentScenarioIndex-th scenario of the selected pool as a
  * Scenario, that is as a std::span< const double >.
  *
  * The function returns a span of the scenarioIndex[currentScenarioIndex]-th 
  * row of the scenarioSet (since we only use the discrete pool approach).
  * 
  * @return A span representing the current scenario
  * @throws std::out_of_range If currentScenarioIndex is out of range
  */
 Scenario get_current_scenario( void ) override;

 /// Function to query the probability weight of the current scenario
 /** When sampling a pool, what are the weights of the drawn scenarios?
  * When using get_scenario_probabilities, we return "the" probability weight
  * inside the pool and not the input probability weight.
  *
  * For scenario pools: We take the input weight and normalize it. Given a
  * scenario with an input_weight, we compute its new_weight by:
  *  "new_weight = input_weight / sum( input_weights_in_the_pool )".
  *
  * This normalization ensures that the weights of all scenarios in the 
  * pool sum to 1.0.
  * 
  * @return The normalized probability of the current scenario
  * @throws std::out_of_range If currentScenarioIndex is out of range
  */
 double get_current_scenario_probability( void ) override;

 /// Move currentScenarioIndex to the next scenario
 /** The function increments currentScenarioIndex by 1 if there is still
  * a scenario left in the pool and returns true. Otherwise, it returns false.
  * 
  * @return true if successfully moved to the next scenario, false if at the end
  */
 bool next_scenario( void ) override;

 /// return the dimension of the scenarios
 /** Every scenario (vector in some Euclidean space R^d) is assumed to have
  * the same dimension. The dimension has been saved in scenarioSize when
  * deserializing the input discrete distribution. */
 ScenarioSize get_scenario_size( void ) override;

/** @} ---------------------------------------------------------------------*/
/*-------------------- SCENARIO POOL MANAGEMENT METHODS --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Scenario Pool Management Methods
 *  @{ */

 /**
  * @brief Randomly select scenarios from the available set
  * 
  * Creates a pool by randomly selecting scenarios with uniform probability.
  * This provides an unbiased sample suitable for simulation and validation
  * of decisions made using a smaller representative set.
  * 
  * @param pool_size Number of scenarios to randomly select
  * @throws std::invalid_argument If pool_size exceeds the available scenarios
  */
 void init_random_pool(ScenarioIndex pool_size) override;

 /**
  * @brief Select the most representative scenarios using optimization
  * 
  * Creates a pool of k scenarios that best approximate the full distribution
  * by minimizing the Wasserstein distance. This produces a high-quality
  * subset for building optimization models that capture the essential
  * characteristics of the uncertainty.
  * 
  * The method can use different algorithms (configurable via set_scenario_reduction_config):
  * - "Dupacova": Fast forward selection with good quality
  * - "BestFit": Better quality through local search
  * - "FirstFit": Balanced speed and quality
  * - "MILP": Optimal selection (slowest)
  * 
  * If no configuration is provided, uses Dupacova algorithm by default.
  * 
  * @param k Number of representative scenarios to select
  * @throws std::invalid_argument If k is invalid (0 or > available scenarios)
  * @throws std::runtime_error If the optimization fails
  * @see set_scenario_reduction_config() to customize the selection algorithm
  */
 void init_representative_pool( ScenarioIndex k ) override;

/** @} ---------------------------------------------------------------------*/
/*----------------- SCENARIO REDUCTION CONFIG METHODS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Scenario Reduction Configuration Methods
 *  @{ */

 /**
  * @brief Get the block configuration for scenario reduction
  * 
  * Returns the BlockConfig part of the scenario reduction configuration,
  * which contains parameters like:
  * - k: Number of scenarios to select
  * - ell: Power parameter for Wasserstein distance (typically 2.0)
  * 
  * The BlockConfig is owned by the DiscreteScenarioSet object and should
  * not be deleted by the caller.
  * 
  * @return Pointer to the BlockConfig, or nullptr if no configuration exists
  */
 BlockConfig* get_scenario_reduction_block_config() const;

 /**
  * @brief Get the solver configuration for scenario reduction
  * 
  * Returns the BlockSolverConfig part of the scenario reduction configuration,
  * which contains parameters like:
  * - algorithm: Scenario reduction method to use (e.g., "Dupacova", "BestFit")
  * - rho: Initial solution parameter (0.0 for random, 1.0 for Dupačová initialization)
  * - shuffle: Whether to shuffle scenarios (for FirstFit algorithm)
  * - Additional solver-specific settings
  * 
  * The BlockSolverConfig is owned by the DiscreteScenarioSet object and should
  * not be deleted by the caller.
  * 
  * @return Pointer to the BlockSolverConfig, or nullptr if no configuration exists
  */
 BlockSolverConfig* get_scenario_reduction_solver_config() const;

 /**
  * @brief Set the scenario reduction configuration
  * 
  * Configures how scenario reduction will be performed when init_representative_pool()
  * is called. Takes ownership of the provided configuration objects.
  * 
  * The provided configurations should contain:
  * 
  * 1. block_config (BlockConfig):
  *    - k: Number of scenarios to select (must be > 0 and <= nbScenarios)
  *    - ell: Power parameter for Wasserstein distance (typically 2.0)
  * 
  * 2. solver_config (BlockSolverConfig):
  *    - algorithm: Scenario reduction method to use, one of:
  *      - "Dupacova" (default, forward selection method)
  *      - "BestFit" (local search with best improvement)
  *      - "FirstFit" (local search with first improvement)
  *      - "MILP" (mixed integer linear programming)
  *    - Other solver-specific parameters
  * 
  * Implementation details:
  * - Stores the configurations internally
  * - Replaces any existing configuration
  * - The object takes ownership of both config pointers
  * 
  * @param block_config BlockConfig containing reduction parameters (k, ell)
  * @param solver_config BlockSolverConfig containing algorithm and solver settings
  */
 void set_scenario_reduction_config(BlockConfig* block_config, BlockSolverConfig* solver_config);
 
 /**
  * @brief Set the scenario reduction configuration with k parameter
  * 
  * Convenience method that sets both the configuration objects and the k parameter.
  * This ensures that k_value is properly set for scenario reduction.
  * 
  * @param block_config BlockConfig containing reduction parameters 
  * @param solver_config BlockSolverConfig containing algorithm and solver settings
  * @param k Number of scenarios to select (must be > 0 and <= nbScenarios)
  */
 void set_scenario_reduction_config(BlockConfig* block_config, BlockSolverConfig* solver_config, ScenarioIndex k);

 /**
  * @brief Set configuration for the DiscreteScenarioSet
  * 
  * This method overrides the base class set_config() to handle scenario
  * reduction configuration. The Configuration object should contain:
  * - k (optional): Number of scenarios to select. If provided and > 0,
  *   init_representative_pool(k) will be called.
  * - BlockConfig (optional): Configuration for the scenario reduction problem
  * - BlockSolverConfig (optional): Configuration for the solver
  * 
  * If k is provided but BlockConfig or BlockSolverConfig are missing,
  * defaults will be generated.
  * 
  * The method takes ownership of any configuration objects extracted from
  * the Configuration.
  * 
  * @param config The Configuration object containing scenario reduction parameters
  */
 void set_config( Configuration* config ) override;

/** @} ---------------------------------------------------------------------*/
/*---------------- METHODS FOR ScenarioReduction FIELDS --------------------*/
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
 
 /// Get the k parameter for scenario reduction
 [[nodiscard]] ScenarioIndex get_k_value() const { return k_value; }
 
 /// Set the k parameter for scenario reduction
 void set_k_value(ScenarioIndex k) {
     if (k == 0) {
         throw std::invalid_argument("k_value must be positive");
     }
     k_value = k;
 }
 
 /// Access an individual scenario value
 [[nodiscard]] double get_scenario_value(ScenarioIndex scenario_idx, ScenarioSize component_idx) const {
     if (scenario_idx >= nbScenarios || component_idx >= scenarioSize) {
         throw std::out_of_range("Index out of range in get_scenario_value");
     }
     return scenarioSet[scenario_idx][component_idx];
 }
 
 /// Get the number of selected scenarios
 [[nodiscard]] size_t get_selected_scenario_count() const {
     return scenarioIndexes.size();
 }
 
 /// Get a specific selected scenario index
 [[nodiscard]] ScenarioIndex get_selected_scenario_index(size_t index) const {
     if (index >= scenarioIndexes.size()) {
         throw std::out_of_range("Index out of range in get_selected_scenario_index");
     }
     return scenarioIndexes[index];
 }

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/** Following members are protected rather than private to allow:
 *  1. Access from test code when validating scenario reduction functionality
 *  2. Access from derived scenario reduction methods that need to manipulate
 *     the scenario sets and selected indexes directly
 */

protected:
 /// Container for Scenario-s
 DiscreteScenarioPool scenarioSet;
 
 /// Indexes of the discrete pool
 std::vector< ScenarioIndex > scenarioIndexes;

 // Use base class deserialize and serialize from the public section
 


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

 /// Pool size
 /** The variable poolSize is initialized when
  * init_representative_pool() or init_random_pool() is used. */
 ScenarioIndex poolSize = 0;

 /// Random generator
 std::mt19937 rng;
 
 /// Flag to check if pool is initialized
 bool is_initialized{false};
 
 /// Compile-time constants
 static constexpr double DEFAULT_EPSILON = 1e-10;
 static constexpr unsigned long DEFAULT_SEED = 1337;
 static constexpr float DEFAULT_ELL_VALUE = 2.0f;
 static constexpr double DEFAULT_RHO_VALUE = 0.0;
 
 /// Number of scenarios to select for scenario reduction
 ScenarioIndex k_value = 0;
 
 /// Pending k value to apply after deserialization
 ScenarioIndex pending_k_value = 0;
 
 /**
  * @brief Configuration for scenario reduction
  * 
  * Stores a pair of configurations for scenario reduction:
  * - First: BlockConfig with parameters like k (number of scenarios) and ell (distance power)
  * - Second: BlockSolverConfig with the algorithm choice and solver settings
  * 
  * The configuration is loaded during deserialization if available
  * in the netCDF file, or can be set programmatically using
  * set_scenario_reduction_config().
  */
 std::pair<BlockConfig*, BlockSolverConfig*> f_scenario_reduction_config = {nullptr, nullptr};

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROBABILITY FIELDS ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Probability weight fields
   * These fields store information about scenario probabilities
 * @{ */

 /// Probabilities of scenarios in the input pool
 /** Vector containing the probability weights of all input scenarios
  * (before any selection). Used to compute normalized probabilities
  * for selected scenarios. */
 std::vector< double > poolProbabilities;

/** @name Fields for the discrete pool
   * When the scenario pool is made of a discrete subset of the input scenarios,
   * it is characterized by a std::vector< ScenarioIndex > and the probability
   * weights can be deduced from the input weights saved in
   * scenarioProbabilities and the sumPoolWeights of the scenarios that
   * belong to the scenario pool.
 * @{ */

 /// holder for the sum of the weights inside the discrete pool. */
 double sumPoolWeights;

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

 /// Empty the scenario pool
 /** Function to clear the internal pool of selected scenario indices.
  * Resets scenarioIndexes, sumPoolWeights, currentScenarioIndex,
  * poolSize, and is_initialized.
  */
 void empty_pool();
 
 /**
  * @brief Extract k parameter with validation
  * 
  * Helper method to safely extract the k parameter (number of scenarios to select)
  * from a BlockConfig, with validation and default value handling.
  * 
  * This method:
  * 1. Extracts the "k" parameter from the BlockConfig
  * 2. Validates that it's positive and doesn't exceed nbScenarios
  * 3. Returns the validated k parameter
  * 
  * Implementation notes:
  * - Throws an exception if k is invalid (not positive or too large)
  * - This is especially important as k controls how many scenarios are selected
  * 
  * @param config The BlockConfig containing the k parameter
  * @return The validated value of k
  * @throws std::runtime_error If the k parameter is invalid (not positive or too large)
  */
 ScenarioIndex get_k_parameter(const BlockConfig* config) const;
 
 /// Determines if scenario reduction should be used based on configuration
 /** This function checks if scenario reduction is configured properly and
  * can be used for the current situation. It evaluates several conditions:
  * 
  * 1. The requested scenario pool size must be greater than 1 (trivial case)
  * 2. A valid scenario reduction configuration must exist
  * 3. The configuration must have valid k and ell parameters
  *
  * This method is called by init_representative_pool() before attempting to use
  * scenario reduction algorithms.
  *
  * @param size The desired size of the reduced scenario pool
  * @return true if scenario reduction can be used, false otherwise
  */
 bool should_use_scenario_reduction(ScenarioIndex size) const;
 
 /// Helper method to apply scenario reduction
 /** This method applies scenario reduction to select a representative subset of scenarios
  * using the configured reduction method. It:
  * 
  * 1. Creates a CapacitatedFacilityLocationBlock to formulate the selection problem
  * 2. Extracts configuration parameters (ell, k, algorithm) from the configuration
  * 3. Configures and runs the ScenarioReductionSolver with the selected algorithm
  * 4. Updates scenarioIndexes with the selected scenario indices
  * 5. Recalculates sumPoolWeights based on the selected scenarios
  * 
  * Implementation notes:
  * - This is called internally by init_representative_pool()
  * - Uses CapacitatedFacilityLocationBlock to formulate the selection problem
  * - Handles different algorithms: Dupacova, BestFit, FirstFit, and MILP
  *
  * @throws std::runtime_error If reduction fails due to configuration or solver issues
  */
 void apply_scenario_reduction();

 /// Validate the k parameter for scenario reduction
 /** Ensures that the k parameter is valid:
  * - Must be greater than 0
  * - Must be less than or equal to the total number of scenarios
  * 
  * @param k The number of scenarios to select
  * @throws std::invalid_argument If k is invalid
  */
 void validate_k_parameter(ScenarioIndex k) const;

 /// Ensure scenario reduction configuration exists
 /** Checks if configuration exists, and creates it with default values if not.
  * This is used by init_representative_pool when no configuration was loaded
  * from deserialization.
  * 
  * @param k The number of scenarios to select
  */
 void ensure_configuration_exists(ScenarioIndex k);

 /// Create CFL problem data structures
 /** Sets up the capacity, fixed cost, and demand vectors for the
  * Capacitated Facility Location problem formulation.
  * 
  * @param n_scenarios The total number of scenarios
  * @return A tuple containing (capacities, fixed_costs, demands)
  */
 std::tuple<CapacitatedFacilityLocationBlock::DVector,
            CapacitatedFacilityLocationBlock::CVector,
            CapacitatedFacilityLocationBlock::DVector>
 create_cfl_problem_data(ScenarioIndex n_scenarios) const;

 /// Compute the transport cost matrix between scenarios
 /** Calculates the distance matrix between all pairs of scenarios using
  * the specified norm (ell parameter).
  * 
  * @param n_scenarios The total number of scenarios
  * @param scenario_size The dimension of each scenario
  * @param ell The power parameter for the distance calculation
  * @return The transport cost matrix
  */
 CapacitatedFacilityLocationBlock::CMatrix
 virtual compute_transport_cost_matrix(ScenarioIndex n_scenarios, 
                              ScenarioSize scenario_size,
                              float ell) const;

 /// Compute distance between two scenarios
 /** Helper method to compute the ell-norm distance between two scenarios.
  * 
  * @param scenario1 First scenario as an Eigen vector
  * @param scenario2 Second scenario as an Eigen vector  
  * @param ell The power parameter for the distance calculation
  * @return The ell-power of the norm distance
  */
 virtual double compute_scenario_distance(const Eigen::VectorXd& scenario1,
                                 const Eigen::VectorXd& scenario2,
                                 float ell) const;

 /// Create and configure the scenario reduction solver
 /** Sets up the ScenarioReductionSolver with the appropriate configuration.
  * 
  * @param cflBlock The configured CapacitatedFacilityLocationBlock
  * @param ell The power parameter for distance calculations
  * @return A configured ScenarioReductionSolver
  */
 Solver*
 create_and_configure_solver(CapacitatedFacilityLocationBlock* cflBlock,
                            float ell) const;

 /// Extract selected scenarios from solver results
 /** Processes the solver results to populate scenarioIndexes with the
  * selected scenario indices.
  * 
  * @param solver The solved ScenarioReductionSolver
  * @param n_scenarios The total number of scenarios
  * @throws std::runtime_error If no scenarios were selected
  */
 void extract_selected_scenarios(const Solver* solver,
                                const CapacitatedFacilityLocationBlock* cflBlock,
                                ScenarioIndex n_scenarios);

 /// Update pool weights after scenario selection
 /** Recalculates sumPoolWeights based on the selected scenarios in
  * scenarioIndexes.
  */
 void update_pool_weights();

 /// Apply pending scenario reduction after deserialization
 /** If a k value was provided during deserialization, this method
  * will call init_representative_pool(k) to apply scenario reduction.
  * This is called at the end of deserialize() to ensure all data
  * is loaded before reduction is attempted.
  */
 void apply_pending_scenario_reduction();


 /// Generate default BlockConfig for CFL
 /** Creates a default BlockConfig suitable for CapacitatedFacilityLocationBlock
  * with the specified k and ell parameters.
  * 
  * @param k Number of scenarios to select
  * @param ell Power parameter for Wasserstein distance (default 2.0)
  * @return A newly allocated BlockConfig (caller owns the pointer)
  */
 BlockConfig* generate_default_cfl_config(ScenarioIndex k, float ell = 2.0) const;

 /// Generate default BlockSolverConfig for ScenarioReductionSolver
 /** Creates a default BlockSolverConfig for the ScenarioReductionSolver
  * with the specified algorithm.
  * 
  * @param algorithm Algorithm name (default "Dupacova")
  * @return A newly allocated BlockSolverConfig (caller owns the pointer)
  */
 BlockSolverConfig* generate_default_solver_config(const std::string& algorithm = "Dupacova") const;

  SMSpp_insert_in_factory_h;

};   // end( class DiscreteScenarioSet )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 } // end( namespace SMSpp_di_unipi_it )

#endif /* __DiscreteScenarioSet */

/*--------------------------------------------------------------------------*/
/*------------------- End file DiscreteScenarioSet.h -----------------------*/
/*--------------------------------------------------------------------------*/