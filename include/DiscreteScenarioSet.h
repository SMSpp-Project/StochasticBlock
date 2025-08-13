/*--------------------------------------------------------------------------*/
/*-------------------- File DiscreteScenarioSet.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class DiscreteScenarioSet, an implementation
 * of ScenarioGenerator for discrete probability distributions stored in netCDF
 * files.
 *
 * DiscreteScenarioSet manages collections of scenarios loaded from netCDF files
 * and provides scenario selection methods including random sampling
 * and Wasserstein distance-based scenario reduction. When configured with
 * a BlockSolverConfig, it formulates the scenario reduction problem as a
 * CapacitatedFacilityLocationBlock instance. The user is responsible for
 * ensuring the chosen Solver is capable of solving this CFL optimization 
 * problem.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Benoît Tran \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Benoît Tran
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
#include "Block.h"           // For BlockConfig
#include "BlockSolverConfig.h" // For BlockSolverConfig
#include "CapacitatedFacilityLocationBlock.h"

#include <Eigen/Dense>  // For Eigen::VectorXd, Eigen::Map
#include <random>       // For std::mt19937
#include <memory>       // For std::unique_ptr

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS DiscreteScenarioSet --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// concrete ScenarioGenerator for discrete probability distributions
/** The DiscreteScenarioSet class is a concrete implementation of ScenarioGenerator
 * that manages discrete probability distributions represented as collections of
 * scenario vectors.
 *
 * ### Data Management
 *
 * Scenarios are loaded from netCDF files and stored internally as a
 * boost::multi_array<double, 2> where each row represents a scenario vector.
 * The class supports both uniform and weighted probability distributions over
 * the scenarios.
 *
 * ### Scenario Selection Methods
 *
 * DiscreteScenarioSet provides two primary methods for selecting scenario subsets:
 *
 * - **Random Selection** (init_random_pool()): Randomly samples \c poolSize scenarios from
 *   the full set, always available without additional configuration
 *
 * - **Scenario Reduction** (init_representative_pool()): Selects a representative
 *   subset that minimizes the Wasserstein distance between the original and reduced
 *   distributions. This requires configuration via set_config() or automatically
 *   loaded during netCDF deserialization
 *
 * ### Configuration and Persistence
 *
 * The class supports configuration through:
 * - Direct parameter setting via set_config() overloads
 * - Configuration objects via set_config()
 * - Persistence through netCDF serialization/deserialization
 *
 * Scenario reduction parameters and solver configurations can be saved to and
 * loaded from netCDF files, enabling reproducible scenario selection.
 *
 * ### Usage Pattern
 *
 * 1. Load scenarios via deserialize() or direct construction
 * 2. Configure scenario reduction (optional) via set_config() overloads
 * 3. Initialize pool via init_random_pool() or init_representative_pool()
 * 4. Access scenarios via get_current_scenario() and iterate via
 *    get_next_scenario()
 *
 * @see ScenarioGenerator for the base interface
 * @see CapacitatedFacilityLocationBlock for the optimization model used
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
 
 /// Type to represent a scenario with its probability
 using ScenarioWithProbability = std::pair<Scenario, double>;

/** @} ---------------------------------------------------------------------*/
/*----------- CONSTRUCTING AND DESTRUCTING DiscreteScenarioSet -------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing DiscreteScenarioSet
 *  @{ */

 DiscreteScenarioSet();

 /// Deserialize the DiscreteScenarioSet from a netCDF group
 /** Read the discrete scenario set from a netCDF::NcGroup. The method expects
  * the following format in the group:
  * 
  * Required dimensions:
  * - "NumberScenarios": Dimension defining the number of scenarios in the pool
  * - "ScenarioSize": Dimension defining the size of each scenario vector
  * 
  * Required variables:
  * - "Scenarios": 2D variable of dimensions [NumberScenarios][ScenarioSize] 
  *   containing the scenario data
  * 
  * Optional variables:
  * - "poolWeights": 1D variable of dimension [NumberScenarios] containing the
  *   probability weights of each scenario. If not present, uniform weights
  *   are assumed. The weights must sum to 1.0 (with tolerance of 1e-6).
  * 
  * Optional subgroup:
  * - "ScenarioReductionConfig": A subgroup containing scenario reduction
  *   configuration, which may include:
  *   - "poolSize": Attribute (scalar) specifying the target pool size for 
  *     scenario reduction
  *   - "ell": Attribute (scalar) specifying the power parameter for the 
  *     ell-Wasserstein distance (default: 2.0)
  *   - "BlockConfig": Subgroup containing serialized BlockConfig for the CFL 
  *     problem
  *   - "BlockSolverConfig": Subgroup containing serialized BlockSolverConfig 
  *     for solving the scenario reduction problem
  * 
  * @param group The netCDF group from which to read the data
  * @throws std::invalid_argument If required dimensions or variables are missing
  * @throws std::runtime_error If data cannot be read from the netCDF file
  * @note This method clears any existing data before loading new data
  */
 void deserialize( const netCDF::NcGroup & group ) override;

 /// Serialize the DiscreteScenarioSet to a netCDF group
 /** Write the discrete scenario set to a netCDF::NcGroup. The method writes
  * the following format to the group:
  * 
  * Created dimensions:
  * - "NumberScenarios": Dimension set to the current number of scenarios
  * - "ScenarioSize": Dimension set to the size of each scenario vector
  * 
  * Created variables:
  * - "Scenarios": 2D variable of dimensions [NumberScenarios][ScenarioSize] 
  *   containing all scenario data
  * - "poolWeights": 1D variable of dimension [NumberScenarios] containing the
  *   probability weights of each scenario
  * 
  * If scenario reduction configuration exists, creates subgroup:
  * - "ScenarioReductionConfig": Subgroup containing:
  *   - "poolSize": Attribute (scalar) with the configured pool size
  *   - "ell": Attribute (scalar) with the ell parameter for distance calculation
  *   - "BlockConfig": Subgroup with serialized BlockConfig if available
  *   - "BlockSolverConfig": Subgroup with serialized BlockSolverConfig if available
  * 
  * @param group The netCDF group to which to write the data
  * @note All scenarios are written, not just the selected pool
  */
 void serialize(netCDF::NcGroup& group) const override;

 virtual ~DiscreteScenarioSet();

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS INHERITED FROM ScenarioGenerator.h ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Functions inherited from ScenarioGenerator.h
 *  @{ */

 /// Set the random seed for reproducible scenario selection
 /** Use this to ensure reproducible results when using random sampling methods.
  * 
  * @param seed The seed value for the random number generator
  */
 void set_seed( unsigned long seed ) override;

 /// Get the current scenario in the iteration
 /** Returns the scenario data at the current position in the selected pool.
  * Use this together with next_scenario() to iterate through all selected scenarios.
  * 
  * @return A read-only view of the current scenario's data
  * @throws std::out_of_range If called before pool initialization or after iteration ends
  * @see next_scenario() to advance to the next scenario
  * @see init_random_pool() or init_representative_pool() to initialize the pool
  */
 Scenario get_current_scenario( void ) const override;

 /// Get the probability weight of the current scenario
 /** Returns the normalized probability of the current scenario within the selected pool.
  * The probabilities of all selected scenarios sum to 1.0.
  * 
  * @return The normalized probability (between 0 and 1)
  * @throws std::runtime_error If pool not initialized
  * @note This is the normalized weight within the pool, not the original weight
  */
 double get_current_scenario_probability( void ) const override;

 /// Move to the next scenario in the iteration
 /** Advances to the next scenario in the selected pool.
  * 
  * @return true if there is a next scenario, false if iteration is complete
  * @see get_current_scenario() to access the scenario data
  */
 bool next_scenario( void ) override;

 /// Get the dimension of scenario vectors
 /** Returns the size of each scenario vector (all scenarios have the same dimension).
  * 
  * @return The number of components in each scenario
  */
 ScenarioSize get_scenario_size( void ) const override;

/** @} ---------------------------------------------------------------------*/
/*-------------------- SCENARIO POOL MANAGEMENT METHODS --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Scenario Pool Management Methods
 *  @{ */

 /// Randomly sample scenarios for Monte Carlo simulation
 /** Creates a random subset of scenarios using weighted sampling.
  * Scenarios with higher weights are more likely to be selected.
  * 
  * @param pool_size Number of scenarios to select (must be ≤ total scenarios)
  * @throws std::invalid_argument If pool_size is 0 or exceeds available scenarios
  * 
  * Example:
  * \code
  * set.init_random_pool(100);  // Select 100 random scenarios
  * while(set.next_scenario()) {
  *     auto scenario = set.get_current_scenario();
  *     // Process scenario...
  * }
  * \endcode
  */
 void init_random_pool(ScenarioIndex pool_size) override;

 /// Select representative scenarios using scenario reduction
 /** Chooses a subset of scenarios that best represents the full distribution.
  * 
  * Two methods are available:
  * - **With solver configuration**: Uses optimization to minimize the Wasserstein
  *   distance between the original and reduced distributions
  * - **Without solver (baseline)**: Selects scenarios with highest weights
  * 
  * Use this when you need:
  * - A smaller set that preserves statistical properties
  * - To reduce computational cost while maintaining solution quality
  * - Scenarios for robust optimization approaches
  * 
  * @param target_pool_size Number of scenarios to select
  * @throws std::invalid_argument If target_pool_size is invalid
  * @throws std::runtime_error If optimization solver fails
  * 
  * Example with optimization:
  * \code
  * // Configure solver for optimal selection
  * set.set_config(block_config, solver_config);
  * set.init_representative_pool(50);  // Find 50 best scenarios
  * \endcode
  * 
  * Example without optimization (baseline):
  * \code
  * set.init_representative_pool(50);  // Select 50 highest-weight scenarios
  * \endcode
  */
 void init_representative_pool( ScenarioIndex target_pool_size ) override;

/** @} ---------------------------------------------------------------------*/
/*----------------- SCENARIO REDUCTION CONFIG METHODS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Scenario Reduction Configuration Methods
 *  @{ */

 /// Check the current scenario reduction configuration
 /** Returns the block configuration used for scenario reduction.
  * 
  * @return The configuration object, or nullptr if not configured
  * @see set_config() to set the configuration
  */
 const BlockConfig* get_block_config() const { return f_block_config; }

 /// Check the solver configuration for scenario reduction
 /** Returns the solver configuration if available.
  * 
  * @return The solver configuration, or nullptr if:
  *         - Not configured yet
  *         - Already used (ownership transferred to solver)
  * @see set_config() to set the configuration
  */
 const BlockSolverConfig* get_solver_config() const { return f_solver_config.get(); }

 /// Configure optimization-based scenario reduction
 /** Sets up the solver and parameters for advanced scenario reduction.
  * When configured, init_representative_pool() will use optimization to find
  * the best representative scenarios.
  * 
  * @param block_config Configuration parameters (will be copied)
  * @param solver_config Solver to use (ownership transferred - do not delete)
  * 
  * Example:
  * \code
  * auto* solver_config = new MILPSolverConfig();
  * set.set_config(block_config, solver_config);
  * // Do NOT delete solver_config - ownership transferred
  * \endcode
  */
 void set_config(BlockConfig* block_config, BlockSolverConfig* solver_config);
 
 /// Configure scenario reduction with explicit pool size
 /** Convenience method to set configuration and pool size together.
  * 
  * @param block_config Configuration parameters (will be copied)
  * @param solver_config Solver configuration (ownership transferred)
  * @param poolSize Number of scenarios to select
  * @throws std::invalid_argument If poolSize is invalid
  */
 void set_config(BlockConfig* block_config, BlockSolverConfig* solver_config, ScenarioIndex poolSize);

 /// Set configuration using a Configuration object
 /** Alternative way to configure scenario reduction using SMS++ Configuration objects.
  * 
  * Supported patterns:
  * - SimpleConfiguration<int>: Just pool size (uses baseline method)
  * - SimpleConfiguration<pair<int,Configuration*>>: Pool size + solver
  * - SimpleConfiguration<pair<Configuration*,Configuration*>>: Full configuration
  * 
  * @param config The configuration object
  * @throws std::invalid_argument If configuration type is not supported
  */
 void set_config( Configuration* config ) override;

/** @} ---------------------------------------------------------------------*/
/*----------------------- GETTERS AND SETTERS ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Getters and Setters
 * Methods to access properties and configuration parameters
 *  @{ */

 /// Get the total number of available scenarios
 /** @return Number of scenarios in the full dataset */
 const ScenarioIndex & get_nbScenarios() const;

 /// Get the dimension of scenario vectors
 /** @return Size of each scenario vector */
 const ScenarioSize & get_scenarioSize() const;
 
 /// Check if a scenario pool has been selected
 /** @return true if init_random_pool() or init_representative_pool() has been called */
 [[nodiscard]] bool is_pool_initialized() const;
 
 /// Get both scenario data and probability together
 /** Convenience method for getting scenario and its weight in one call.
  * @return Pair of (scenario data, normalized probability) */
 [[nodiscard]] ScenarioWithProbability get_current_scenario_with_prob() const;
 
 /// Access a specific value from any scenario
 /** Direct access to individual scenario components.
  * @param scenario_idx Which scenario (0 to nbScenarios-1)
  * @param component_idx Which component of that scenario (0 to scenarioSize-1)
  * @return The requested value
  * @throws std::out_of_range If indices are invalid */
 [[nodiscard]] double get_scenario_value(ScenarioIndex scenario_idx, ScenarioSize component_idx) const {
     if (scenario_idx >= nbScenarios || component_idx >= scenarioSize) {
         throw std::out_of_range("Index out of range in get_scenario_value");
     }
     return scenarioSet[scenario_idx][component_idx];
 }
 
 /// Get all indices of selected scenarios
 /** Returns which scenarios were selected by init_random_pool() or init_representative_pool().
  * @return Vector of scenario indices in the selected pool
  * @throws std::runtime_error If pool not initialized */
 const std::vector<ScenarioIndex>& get_selected_scenarios() const;
 
 /// Get the index of a specific selected scenario
 /** @param index Position in the selected pool (0 to poolSize-1)
  * @return The original index of that scenario in the full dataset
  * @throws std::out_of_range If index >= pool size */
 [[nodiscard]] ScenarioIndex get_selected_scenario_index(size_t index) const {
     if (index >= scenarioIndexes.size()) {
         throw std::out_of_range("Index out of range in get_selected_scenario_index");
     }
     return scenarioIndexes[index];
 }

 /// Get the configured pool size
 /** @return Number of scenarios that will be selected */
 [[nodiscard]] ScenarioIndex get_poolSize() const { return poolSize; }
 
 /// Get the distance power parameter
 /** @return The ell parameter for ell-Wasserstein distance (default: 2.0) */
 [[nodiscard]] float get_ell() const { return ell; }
 
 /// Set the distance power parameter
 /** Controls how distances are calculated in scenario reduction.
  * Common values: 1.0 (linear), 2.0 (quadratic, default)
  * @param ell_value Power parameter (must be > 0)
  * @throws std::invalid_argument If ell_value ≤ 0 */
 void set_ell(float ell_value) {
     if (ell_value <= 0) {
         throw std::invalid_argument("ell must be positive");
     }
     ell = ell_value;
 }

/** @} ---------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

private:

 /// Container for Scenario-s
 DiscreteScenarioPool scenarioSet;
 
 /// Indexes of the discrete pool
 std::vector< ScenarioIndex > scenarioIndexes;

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
 /** The variable \c poolSize is initialized when
  * \c init_representative_pool() or \c init_random_pool() is used. */
 ScenarioIndex poolSize = 0;

 /// Random generator
 std::mt19937 rng;
 
 /// Flag to check if pool is initialized
 bool is_initialized{false};
 
 /// Power parameter for the ell-Wasserstein distance in scenario reduction
 /** The \c ell parameter determines the power used in the ell-Wasserstein distance
  *  calculation for scenario reduction. This is the power to which the ground metric
  *  (typically Euclidean distance) is raised.
  *  
  *  In the transportation cost computation:
  *  - transport_cost[i][j] = ||scenario_i - scenario_j||^ell
  *  
  *  Common values:
  *  - 1.0: Linear transportation costs (1-Wasserstein distance)
  *  - 2.0: Quadratic transportation costs (2-Wasserstein distance) [DEFAULT]
  *  
  *  The default value of 2.0 is chosen because:
  *  - The 2-Wasserstein distance has favorable theoretical properties
  *  - It provides a good balance between computational efficiency and accuracy
  *  - It is the most commonly used value in scenario reduction literature
  *  
  *  This value can be changed via \c set_ell() or during deserialization from netCDF. */
 float ell = 2.0f;  // Default: 2-Wasserstein distance
 
 /// Block configuration for scenario reduction
 /** Stores the BlockConfig used for creating the CapacitatedFacilityLocationBlock
  * during scenario reduction. */
 BlockConfig* f_block_config = nullptr;
 
 /// Solver configuration for scenario reduction
 /** Stores the BlockSolverConfig used for solving the CFL problem.
  * 
  * This is a unique_ptr that manages the lifetime of the BlockSolverConfig.
  * When set_config() is called, ownership of the provided BlockSolverConfig
  * is transferred to this unique_ptr.
  * 
  * IMPORTANT: During scenario reduction (in create_and_configure_solver()),
  * when apply() is called on the BlockSolverConfig, ownership is transferred
  * to the solver via release(). After this point, f_solver_config becomes
  * nullptr and get_solver_config() will return nullptr.
  * 
  * The field is mutable because ownership transfer happens in the const method
  * create_and_configure_solver() during scenario reduction. */
 mutable std::unique_ptr<BlockSolverConfig> f_solver_config;

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROBABILITY FIELDS ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Probability weight fields
   * These fields store information about scenario probabilities
 * @{ */

 /// Weights of scenarios in the input pool
 /** Vector containing the weights of all input scenarios
  * (before any selection). Used to compute normalized probabilities
  * for selected scenarios. */
 std::vector< double > poolWeights;

/** @name Pool weight management fields
   * These fields manage the weights and probabilities of selected scenarios
 * @{ */

 /// Sum of the weights of scenarios in the discrete pool
 /** Holds the sum of weights for all selected scenarios in \c scenarioIndexes.
  * Used as denominator when computing normalized probabilities. */
 double sumPoolWeights;
 
 /// Normalized weights of scenarios in the current pool
 /** Vector containing the normalized weights of scenarios that have been
  * selected for the current pool (via init_random_pool or init_representative_pool).
  * These weights are computed as: <tt>poolWeights[scenario_i] / sumPoolWeights</tt>.
  * The size of this vector equals the number of selected scenarios (\c poolSize).
  * This container is only populated after pool initialization. */
 std::vector< double > normalizedPoolWeights;

/** @} ---------------------------------------------------------------------*/
/*-------------------- HELPER METHODS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name helper methods of the class
 * Miscellaneous functions
 * @{ */

 /// Empty the scenario pool
 /** Function to clear the internal pool of selected scenario indices.
  * Resets \c scenarioIndexes, \c sumPoolWeights, \c currentScenarioIndex,
  * \c poolSize, and \c is_initialized.
  */
 void empty_pool();


 /// Compute the transport cost matrix between scenarios
 /** Calculates the distance matrix between all pairs of scenarios using
  * the specified norm (ell parameter).
  * 
  * @param n_scenarios The total number of scenarios
  * @param scenario_size The dimension of each scenario
  * @param ell The power parameter for the distance calculation
  * @return The transport cost matrix
  */
 virtual CapacitatedFacilityLocationBlock::CMatrix
 compute_transport_cost_matrix(ScenarioIndex n_scenarios, 
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

 /// Update pool weights after scenario selection
 /** Computes \c sumPoolWeights and populates \c normalizedPoolWeights based on 
  * the selected scenarios in \c scenarioIndexes after either \c init_random_pool() 
  * or \c init_representative_pool() have been used.
  * 
  * This method:
  * 1. Calculates the sum of weights for all selected scenarios (\c sumPoolWeights)
  * 2. Populates the \c normalizedPoolWeights container with normalized probabilities
  *    computed as: <tt>poolWeights[scenario_i] / sumPoolWeights</tt>
  * 3. Falls back to uniform distribution if \c sumPoolWeights is zero.
  */
 void update_pool_weights();
 
 /// Apply baseline selection method
 /** Selects the top scenarios based on their probability weights.
  * This is the fallback method when no BlockSolverConfig is provided.
  * Scenarios with higher weights are selected first.
  * 
  * @param target_size Number of scenarios to select
  */
 void apply_baseline_selection(ScenarioIndex target_size);

  SMSpp_insert_in_factory_h;

};   // end( class DiscreteScenarioSet )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 } // end( namespace SMSpp_di_unipi_it )

#endif /* __DiscreteScenarioSet */

/*--------------------------------------------------------------------------*/
/*------------------- End file DiscreteScenarioSet.h -----------------------*/
/*--------------------------------------------------------------------------*/