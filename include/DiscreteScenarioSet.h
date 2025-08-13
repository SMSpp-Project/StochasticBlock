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
 /** Extends ScenarioGenerator::deserialize( netCDF::NcGroup ) to the
  * specific format of a DiscreteScenarioSet. Besides what is managed by the
  * serialize() method of the base ScenarioGenerator class, the group should
  * contain the following:
  *
  * - the dimension "NumberScenarios" containing the number of scenarios in
  *   the discrete probability distribution (must be positive)
  *
  * - the dimension "ScenarioSize" containing the dimension of the Euclidean
  *   space (R^d) representing a single scenario (must be positive)
  *
  * - the variable "Scenarios", of type double and indexed over both the
  *   dimensions "NumberScenarios" and "ScenarioSize"; the entry ( i , j )
  *   is assumed to contain the j-th component of the i-th scenario
  *
  * The dimensions "NumberScenarios" and "ScenarioSize" and the variable
  * "Scenarios" are mandatory. However, the optional
  *
  * - variable "poolWeights", of type double and indexed over the
  *   dimension "NumberScenarios"; the i-th entry of the variable is assumed
  *   to contain the weight of the i-th scenario
  *
  * may also be present to represent the probabilities of the scenarios (if
  * not present, uniform probabilities 1.0/NumberScenarios are assumed). If 
  * provided, the variable must have exactly NumberScenarios elements and
  * probabilities must sum to approximately 1.0 (within 1e-6 tolerance). Also,
  * the optional
  *
  * - group "ScenarioReductionConfig" can be present; if so, it may contain:
  *
  *   = variable "poolSize", of type int, containing the number of scenarios to select
  *     for scenario reduction (takes precedence over any \c SimpleConfiguration<int>
  *     in BlockConfig's \c extra_Configuration if both exist)
  *
  *   = variable "ell", of type float, containing the power parameter for
  *     Wasserstein distance calculation (default 2.0)
  *
  *   = group "BlockConfig" containing the configuration for the
  *     CapacitatedFacilityLocationBlock used in scenario reduction. If this
  *     contains a \c SimpleConfiguration<int> in its \c extra_Configuration and no
  *     "poolSize" variable exists, the integer value is used as \c poolSize. If BlockConfig
  *     is not provided but \c poolSize is, a default BlockConfig is generated.
  *
  *   = group "BlockSolverConfig" containing the configuration for the solver
  *     used with the CapacitatedFacilityLocationBlock. If not provided, a
  *     default solver configuration is generated.
  *
  * If "ScenarioReductionConfig" is present, \c poolSize must be provided either as a
  * variable or within BlockConfig's \c extra_Configuration, otherwise an error
  * is thrown.
  *
  * Note: This method clears any existing data and configuration before
  * deserializing. 
  */
 void deserialize( const netCDF::NcGroup & group ) override;

 /**
  * @brief Serialize the object to a netCDF group
  * 
  * Extends the base class serialization to include scenario reduction configuration.
  * Writes the DiscreteScenarioSet data in the format expected by deserialize().
  * 
  * This method:
  * 1. Calls the base class serialize (from ScenarioGenerator)
  * 2. Writes mandatory data:
  *    - dimension "NumberScenarios" with the number of scenarios
  *    - dimension "ScenarioSize" with the scenario dimension
  *    - variable "Scenarios" (double array) with all scenario data
  *    - variable "poolWeights" (double array) with scenario weights
  * 3. If scenario reduction configuration exists, creates "ScenarioReductionConfig" group
  * 4. Within "ScenarioReductionConfig", writes:
  *    - variable "poolSize" (int) with the number of scenarios to select
  *    - variable "ell" (float) with the Wasserstein distance power parameter
  * 5. Creates subgroups "BlockConfig" and "BlockSolverConfig" within
  *    "ScenarioReductionConfig" for the respective configurations
  * 
  * Note that scenario reduction configuration is only serialized if it exists
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

 /// Setting the seed of the pseudo-random number generator
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
 Scenario get_current_scenario( void ) const override;

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
 double get_current_scenario_probability( void ) const override;

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
 ScenarioSize get_scenario_size( void ) const override;

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
  * @brief Select the most representative scenarios
  * 
  * Creates a pool of scenarios from the full distribution. The selection
  * method depends on whether a BlockSolverConfig is provided:
  * 
  * 1. If BlockSolverConfig is provided:
  *    - Creates a CapacitatedFacilityLocationBlock to model scenario selection
  *    - Applies BlockConfig if provided (or generates minimal default)
  *    - Applies the configured solver (e.g., "Dupacova", "BestFit", "MILP")
  *    - Solves to minimize Wasserstein distance between original and reduced sets
  * 
  * 2. If no BlockSolverConfig is provided:
  *    - Uses baseline method (selects top scenarios by probability weight)
  * 
  * The baseline method sorts scenarios by their probability weights in descending
  * order and selects the scenarios with highest weights.
  * 
  * @param target_pool_size Number of representative scenarios to select
  * @throws std::invalid_argument If target_pool_size is invalid (0 or > available scenarios)
  * @throws std::runtime_error If the solver fails (when using configured solver)
  * @see set_config() to provide BlockConfig and BlockSolverConfig
  */
 void init_representative_pool( ScenarioIndex target_pool_size ) override;

/** @} ---------------------------------------------------------------------*/
/*----------------- SCENARIO REDUCTION CONFIG METHODS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Scenario Reduction Configuration Methods
 *  @{ */

 /**
  * @brief Get the block configuration for scenario reduction
  * 
  * @return Const pointer to the BlockConfig, or nullptr if no configuration exists
  */
 const BlockConfig* get_block_config() const { return f_block_config; }

 /**
  * @brief Get the solver configuration for scenario reduction
  * 
  * @return Const pointer to the BlockSolverConfig, or nullptr in these cases:
  *         - No configuration has been set yet
  *         - Ownership was already transferred to a solver during scenario reduction
  *           (after init_representative_pool() has been called with a configured solver)
  * 
  * @note After successful scenario reduction using a solver, this will return nullptr
  *       because the BlockSolverConfig ownership is transferred to the solver.
  */
 const BlockSolverConfig* get_solver_config() const { return f_solver_config.get(); }

 /**
  * @brief Set the scenario reduction configuration
  * 
  * Configures how scenario reduction will be performed when init_representative_pool()
  * is called.
  * 
  * Memory management:
  * - block_config: This method creates a clone. You retain ownership of the original.
  * - solver_config: This method takes ownership. Do NOT delete it after this call.
  * 
  * @param block_config BlockConfig containing reduction parameters like \c poolSize (will be cloned)
  * @param solver_config BlockSolverConfig containing solver settings (ownership transferred)
  */
 void set_config(BlockConfig* block_config, BlockSolverConfig* solver_config);
 
 /**
  * @brief Set the scenario reduction configuration with \c poolSize parameter (convenience overload)
  * 
  * Convenience method that sets both the configuration objects and the \c poolSize parameter.
  * This ensures that \c poolSize is properly set for scenario reduction.
  * 
  * @param block_config BlockConfig containing reduction parameters 
  * @param solver_config BlockSolverConfig containing the Solver to use for CFL optimization
  * @param poolSize Number of scenarios to select (must be > 0 and <= nbScenarios)
  */
 void set_config(BlockConfig* block_config, BlockSolverConfig* solver_config, ScenarioIndex poolSize);

 /**
  * @brief Set configuration for the DiscreteScenarioSet
  * 
  * Supports multiple configuration patterns:
  * - \c SimpleConfiguration<int>: \c poolSize only (baseline method)
  * - \c SimpleConfiguration<pair<int,Configuration*>>: \c poolSize + solver config
  * - SimpleConfiguration<pair<Configuration*, Configuration*>>: full block + solver config
  * 
  * @param config The Configuration object containing scenario reduction parameters
  */
 void set_config( Configuration* config ) override;

/** @} ---------------------------------------------------------------------*/
/*------------------- GETTERS AND SETTERS METHODS --------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Getters and Setters for DiscreteScenarioSet properties
 *  @{ */

 // Basic Properties Getters
 /// get a reference to nbScenarios
 const ScenarioIndex & get_nbScenarios() const;

 /// get a reference to scenarioSize
 const ScenarioSize & get_scenarioSize() const;
 
 /// Check if the scenario pool has been initialized
 [[nodiscard]] bool is_pool_initialized() const;
 
 // Scenario Access Methods
 /// Get current scenario with its probability as a pair
 [[nodiscard]] ScenarioWithProbability get_current_scenario_with_prob() const;
 
 /// Access an individual scenario value
 [[nodiscard]] double get_scenario_value(ScenarioIndex scenario_idx, ScenarioSize component_idx) const {
     if (scenario_idx >= nbScenarios || component_idx >= scenarioSize) {
         throw std::out_of_range("Index out of range in get_scenario_value");
     }
     return scenarioSet[scenario_idx][component_idx];
 }
 
 /// Get a specific selected scenario index
 [[nodiscard]] ScenarioIndex get_selected_scenario_index(size_t index) const {
     if (index >= scenarioIndexes.size()) {
         throw std::out_of_range("Index out of range in get_selected_scenario_index");
     }
     return scenarioIndexes[index];
 }

 /// Get the \c poolSize parameter for scenario reduction
 [[nodiscard]] ScenarioIndex get_poolSize() const { return poolSize; }
 
 /// Set the \c poolSize parameter for scenario reduction
 void set_poolSize(ScenarioIndex poolSize);
 
 /// Get the ell parameter for Wasserstein distance calculation
 [[nodiscard]] float get_ell() const { return ell; }
 
 /// Set the ell parameter for Wasserstein distance calculation
 void set_ell(float ell_value) {
     if (ell_value <= 0) {
         throw std::invalid_argument("ell must be positive");
     }
     ell = ell_value;
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
 /** The variable \c poolSize is initialized when
  * \c init_representative_pool() or \c init_random_pool() is used. */
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
 
 /// Power parameter for Wasserstein distance calculation in scenario reduction
 /** The ell parameter determines the power of the norm used when computing
  *  distances between scenarios. Default value is 2.0 (Euclidean distance).
  *  This can be overridden during deserialization from the netCDF file. */
 float ell = DEFAULT_ELL_VALUE;
 
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

/** @name Fields for the discrete pool
   * When the scenario pool is made of a discrete subset of the input scenarios,
   * it is characterized by a std::vector< ScenarioIndex > and the probability
   * weights can be deduced from the input weights saved in
   * \c poolWeights and the \c sumPoolWeights of the scenarios that
   * belong to the scenario pool.
 * @{ */

 /// holder for the sum of the weights inside the discrete pool. */
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


 /// Extract selected scenarios from CFL block results
 /** Reads the y variables from the CFL block to populate \c scenarioIndexes with the
  * selected scenario indices. The solver solution must have been written to the
  * block variables before calling this function.
  * 
  * @param cflBlock The CapacitatedFacilityLocationBlock containing the solution
  * @param n_scenarios The total number of scenarios
  * @throws std::runtime_error If no scenarios were selected
  */
 void get_selected_scenarios_from_block(const CapacitatedFacilityLocationBlock* cflBlock,
                                       ScenarioIndex n_scenarios);

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


//  /// Generate default BlockSolverConfig for ScenarioReductionSolver
//  /** Creates a default BlockSolverConfig for the ScenarioReductionSolver
//   * with the specified algorithm.
//   * 
//   * @param algorithm Algorithm name (default "Dupacova")
//   * @return A newly allocated BlockSolverConfig (caller owns the pointer)
//   */
//  BlockSolverConfig* generate_default_solver_config(const std::string& algorithm = "Dupacova") const;
 
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