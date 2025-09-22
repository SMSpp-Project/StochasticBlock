/*--------------------------------------------------------------------------*/
/*-------------------- File DiscreteScenarioSet.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class DiscreteScenarioSet, an
 * implementation of ScenarioGenerator for discrete probability
 * distributions. Thus, a DiscreteScenarioSet is basically just a collection
 * of vectors, one for each scenario. The main way in which these are loaded
 * is from netCDF files.
 *
 * DiscreteScenarioSet manages finite collections of scenarios and provides
 * scenario selection methods, including random sampling and Wasserstein
 * distance-based scenario reduction. When configured with a
 * BlockSolverConfig, it formulates the scenario reduction problem as a
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
#include "Block.h"              // for BlockConfig
#include "BlockSolverConfig.h"  // for BlockSolverConfig

#include <Eigen/Dense>  // for Eigen::VectorXd, Eigen::Map
#include <random>       // for std::mt19937

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
/** The DiscreteScenarioSet class is a concrete implementation of
 * ScenarioGenerator that manages discrete probability distributions
 * represented as collections of scenario vectors.
 *
 * ### Data Management
 *
 * Scenarios are loaded from netCDF files and stored internally as a
 * boost::multi_array< double , 2 > where each row represents a scenario
 * vector. The class supports both uniform and weighted probability
 * distributions over the scenarios.
 *
 * ### Scenario Selection Methods
 *
 * DiscreteScenarioSet provides two primary methods for selecting scenario
 * subsets:
 *
 * - **Random Selection** (init_random_pool()): randomly samples \c poolSize
 *   scenarios from the full set, always available without additional
 *   configuration.
 *
 * - **Scenario Reduction** (init_representative_pool()): selects a
 *   representative subset that minimizes the Wasserstein distance between
 *   the original and reduced distributions. This requires an appropriate
 *   Configuration, that can be passed either via set_config() or loaded
 *   during netCDF deserialization
 *
 * ### Configuration and Persistence
 *
 * The class supports configuration through:
 * - Direct parameter setting via set_config()
 * - Configuration objects via set_config()
 * - Persistence through netCDF serialization/deserialization
 *
 * Scenario reduction parameters and solver configurations can be saved to
 * and loaded from netCDF files, enabling reproducible scenario selection.
 *
 * ### Usage Pattern
 *
 * 1. Load scenarios via deserialize() or direct construction
 * 2. Configure scenario reduction (optional) via set_config()
 * 3. Initialize pool via init_random_pool() or init_representative_pool()
 * 4. Access scenarios via get_current_scenario() and iterate via
 *    next_scenario()
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

 /// Type for cost matrix (copied from CapacitatedFacilityLocationBlock)
 /** This type definition matches CapacitatedFacilityLocationBlock::CMatrix
  * to avoid needing to include the full header. */
 
 using CMatrix = boost::multi_array< double , 2 >;

/** @} ---------------------------------------------------------------------*/
/*----------- CONSTRUCTING AND DESTRUCTING DiscreteScenarioSet -------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing DiscreteScenarioSet
 *  @{ */

 DiscreteScenarioSet();

 /// Deserialize the DiscreteScenarioSet from a netCDF group
 /** Read the discrete scenario set from a netCDF::NcGroup. The method
  * expects the following format in the group:
  * 
  * Required dimensions:
  *
  * - "NumberScenarios": Dimension defining the number of scenarios in the
  *   pool
  *
  * - "ScenarioSize": Dimension defining the size of each scenario vector
  * 
  * Required variables:
  *
  * - "Scenarios": 2D variable of dimensions [NumberScenarios][ScenarioSize] 
  *   containing the scenario data
  * 
  * Optional variables:
  *
  * - "poolWeights": 1D variable of dimension [NumberScenarios] containing the
  *   probability weights of each scenario. If not present, uniform weights
  *   are assumed. The weights must sum to 1.0 (with tolerance of 1e-6).
  *
  * @param group The netCDF group from which to read the data
  * @throws std::invalid_argument If required dimensions or variables are
  *              missing
  * @throws std::runtime_error If data cannot be read from the netCDF file
  * @note This method clears any existing data before loading new data.
  * @note ScenarioReductionConfig is NOT loaded from the netCDF file. If
  *       scenario reduction functionality is needed after deserialization,
  *       the configuration must be set separately using set_config(). */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
 /// Serialize the DiscreteScenarioSet to a netCDF group
 /** Write the discrete scenario set to a netCDF::NcGroup with the format
  * described by deserialize().
  *
  * @param group The netCDF group to which to write the data
  * @note All scenarios are written, not just the selected pool
  * @note ScenarioReductionConfig is NOT serialized. If scenario reduction
  *       configuration is needed after deserialization, it must be set
  *       separately using set_config() */

 void serialize( netCDF::NcGroup& group ) const override;

/*--------------------------------------------------------------------------*/

 virtual ~DiscreteScenarioSet();

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS INHERITED FROM ScenarioGenerator.h ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Functions inherited from ScenarioGenerator.h
 *  @{ */

 /// Set the random seed for reproducible scenario selection
 /** Use this to ensure reproducible results when using random sampling
  * methods.
  * 
  * @param seed The seed value for the random number generator */

 void set_seed( unsigned long seed ) override { rng.seed( seed ); }

 /// Get the current scenario in the iteration
 /** Returns the scenario data at the current position in the selected pool.
  * Use this together with next_scenario() to iterate through all selected
  * scenarios.
  * 
  * @return A read-only view of the current scenario's data
  * @throws std::out_of_range If called before pool initialization or after
  *         iteration ends
  * @see next_scenario() to advance to the next scenario
  * @see init_random_pool() or init_representative_pool() to initialize the
  *      pool */

 [[nodiscard]] inline Scenario get_current_scenario( void ) const override {
   const auto index = scenarioIndexes[ currentScenarioIndex ];
   return Scenario( & scenarioSet[ index ][ 0 ] , scenarioSize );
   }

/*--------------------------------------------------------------------------*/
 /// Get the probability weight of the current scenario
 /** Returns the normalized probability of the current scenario within the
  * selected pool. The probabilities of all selected scenarios sum to 1.0.
  * 
  * @return The normalized probability (between 0 and 1)
  * @throws std::runtime_error If pool not initialized
  * @note This is the normalized weight within the pool, not the original
  *       weight */

 [[nodiscard]] double get_current_scenario_probability( void ) const
  override { return( poolWeights[ currentScenarioIndex ] ); }

/*--------------------------------------------------------------------------*/
 /// Move to the next scenario in the iteration
 /** Advances to the next scenario in the selected pool.
  * 
  * @return true if there is a next scenario, false if iteration is complete
  * @see get_current_scenario() to access the scenario data */

 bool next_scenario( void ) override;

/*--------------------------------------------------------------------------*/
 /// Get the dimension of scenario vectors
 /** Returns the size of each scenario vector (all scenarios have the same
  * dimension).
  * 
  * @return The number of components in each scenario */

 [[nodiscard]] ScenarioSize get_scenario_size( void ) const override {
  return( scenarioSize );
  }

/** @} ---------------------------------------------------------------------*/
/*-------------------- SCENARIO POOL MANAGEMENT METHODS --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Scenario Pool Management Methods
 *  @{ */

 /// Randomly sample scenarios for Monte Carlo simulation
 /** Creates a random subset of scenarios using weighted sampling.
  * Scenarios with higher weights are more likely to be selected.
  * 
  * @param pool_size Number of scenarios to select, which must be less than
  *        or equal to the total number of scenarios
  * @throws std::invalid_argument If pool_size is 0 or exceeds available
  *         scenarios
  * 
  * Example:
  * \code
  * set.init_random_pool( 100 );  // Select 100 random scenarios
  * while( set.next_scenario() ) {
  *  auto scenario = set.get_current_scenario();
  *  // Process scenario...
  *  }
  * \endcode
  */

 void init_random_pool(ScenarioIndex pool_size) override;

/*--------------------------------------------------------------------------*/
 /// Select representative scenarios using scenario reduction
 /** Chooses a subset of scenarios that best represents the full distribution.
  * 
  * Two methods are available:
  *
  * - **With solver configuration**: uses optimization to minimize the
  *   Wasserstein distance between the original and reduced distributions

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

 [[nodiscard]] const BlockConfig * get_block_config( void ) const {
  return( f_block_config );
  }

/*--------------------------------------------------------------------------*/
 /// Check the solver configuration for scenario reduction
 /** Returns the solver configuration if set.
  * 
  * @return The solver configuration, or nullptr if not configured
  * @see set_config() to set the configuration */
 
 [[nodiscard]] const BlockSolverConfig * get_solver_config( void ) const {
  return( f_solver_config.get() );
  }

/*--------------------------------------------------------------------------*/
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

 void set_config( BlockConfig * block_config ,
		  BlockSolverConfig * solver_config );

/*--------------------------------------------------------------------------*/
 /// Configure scenario reduction with explicit pool size
 /** Convenience method to set configuration and pool size together.
  * 
  * @param block_config Configuration parameters (will be copied)
  * @param solver_config Solver configuration (ownership transferred)
  * @param poolSize Number of scenarios to select
  * @throws std::invalid_argument If poolSize is invalid
  */

 void set_config( BlockConfig * block_config ,
		  BlockSolverConfig * solver_config ,
		  ScenarioIndex poolSize );

/*--------------------------------------------------------------------------*/
 /// Set configuration using a Configuration object
 /** Alternative way to configure scenario reduction using SMS++
  * Configuration objects.
  * 
  * Supported patterns:
  *
  * - SimpleConfiguration< int >: Just pool size (uses baseline method)
  *
  * - SimpleConfiguration<pair< int , Configuration* > >: Pool size + solver
  *
  * - SimpleConfiguration< pair< Configuration* , Configuration* > >:
  *   Full configuration
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

 [[nodiscard]] const ScenarioIndex & get_nbScenarios( void ) const {
  return( nbScenarios );
  }

/*--------------------------------------------------------------------------*/
 /// Get the dimension of scenario vectors
 /** @return Size of each scenario vector */

 [[nodiscard]] const ScenarioSize & get_scenarioSize( void ) const {
  return( scenarioSize );
  }
 
/*--------------------------------------------------------------------------*/
 /// Check if a scenario pool has been selected
 /** @return true if init_random_pool() or init_representative_pool() has been called */

 [[nodiscard]] bool is_pool_initialized( void ) const override {
  return( is_initialized );
  }

/*--------------------------------------------------------------------------*/
 /// Reset the pool iteration to the beginning
 /** Resets the current scenario index to 0, allowing iteration through the
  * same pool from the beginning without re-initializing.
  *
  * @throws std::runtime_error If no pool has been initialized */

 void reset_pool( void ) override {
  if( ! is_initialized )
   throw( std::runtime_error( "reset_pool: pool not initialized" ) );

  currentScenarioIndex = 0;
  }

/*--------------------------------------------------------------------------*/
 /// Access a specific value from any scenario
 /** Direct access to individual scenario components.
  * @param scenario_idx Which scenario (0 to nbScenarios-1)
  * @param component_idx Which component of that scenario (0 to
  *        scenarioSize-1)
  * @return The requested value
  * @throws std::out_of_range If indices are invalid */

 [[nodiscard]] double get_scenario_value( ScenarioIndex scenario_idx ,
					  ScenarioSize component_idx ) const
  {
   if( scenario_idx >= nbScenarios || component_idx >= scenarioSize )
    throw( std::out_of_range( "get_scenario_value:: index out of range" ) );

   return( scenarioSet[ scenario_idx ][ component_idx ] );
   }

/*--------------------------------------------------------------------------*/
 /// Get all indices of selected scenarios
 /** Returns which scenarios were selected by init_random_pool() or
  * init_representative_pool().
  * @return Read-only view of scenario indices in the selected pool
  * @throws std::runtime_error If pool not initialized */

 [[nodiscard]] std::span< const ScenarioIndex > get_selected_scenarios()
  const;

/*--------------------------------------------------------------------------*/
 /// Get the configured pool size
 /** @return Number of scenarios that will be selected */

 [[nodiscard]] ScenarioIndex get_poolSize( void ) const {
  return( poolSize );
  }
 
/*--------------------------------------------------------------------------*/
 /// Get the distance power parameter
 /** @return The ell parameter for ell-Wasserstein distance (default: 2.0) */

 [[nodiscard]] float get_ell( void ) const { return( ell ); }
 
/*--------------------------------------------------------------------------*/
 /// Get all scenario weights from the input set
 /** Returns the weights of all scenarios in the input set (not just
  * selected ones).
  * @return Read-only view of all scenario weights from the input set */

 [[nodiscard]] inline std::span<const double> get_set_weights( void ) const { 
  return( setWeights );
  }

/*--------------------------------------------------------------------------*/
 /// Get weights of selected scenarios in the pool
 /** Returns the weights of scenarios in the current pool after selection.
  * These weights sum to 1.0 and correspond to the selected scenarios.
  * For init_random_pool: simple renormalization of original weights
  * For init_representative_pool: aggregated weights from assignments
  * @return Read-only view of pool weights
  * @throws std::runtime_error If pool not initialized */

 [[nodiscard]] inline std::span<const double> get_pool_weights( void )
  const {
  if( ! is_initialized )
   throw( std::runtime_error( "get_pool_weights:: pool not initialized" ) );
  return( poolWeights );
  }
 
/*--------------------------------------------------------------------------*/
 /// Set the distance power parameter
 /** Controls how distances are calculated in scenario reduction.
  * Common values: 1.0 (linear), 2.0 (quadratic, default)
  * @param ell_value Power parameter (must be > 0)
  * @throws std::invalid_argument If ell_value ≤ 0 */

 void set_ell( float ell_value ) {
  if( ell_value <= 0 )
   throw( std::invalid_argument( "set_ell: ell must be positive" ) );

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
 /** The \c ell parameter determines the power used in the ell-Wasserstein
  * distance calculation for scenario reduction. This is the power to which
  * the ground metric (default Euclidean distance) is raised.
  *  
  *  In the transportation cost computation:
  *  - transport_cost[i][j] = || scenario_i - scenario_j ||^ell
  *  
  *  Common values:
  *  - 1.0: Linear transportation costs (1-Wasserstein distance)
  *  - 2.0: Quadratic transportation costs (2-Wasserstein distance) [DEFAULT]
  *  
  *  This value can be changed via \c set_ell() or during deserialization
  * from netCDF. */
 float ell = 2.0f;  // Default: 2-Wasserstein distance

 /// Block configuration for scenario reduction
 /** Stores the BlockConfig used for creating the
  * CapacitatedFacilityLocationBlock during scenario reduction. */

 BlockConfig * f_block_config = nullptr;

 /// Solver configuration for scenario reduction
 /** Stores the BlockSolverConfig used for solving the CFL problem.
  * 
  * When init_representative_pool() is called, this configuration is
  * cloned and applied to create a new solver for the
  * CapacitatedFacilityLocationBlock. The original configuration
  * remains available for reuse in subsequent calls. */

 std::unique_ptr< BlockSolverConfig > f_solver_config;

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROBABILITY FIELDS ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Probability weight fields
 * These fields store information about scenario probabilities
 * @{ */

 /// Weights of scenarios in the input set
 /** Vector containing the weights of all scenarios in the input set
  * (before any selection). These are the weights loaded from deserialization.
  * Used as input for computing pool weights after selection. */
 std::vector< double > setWeights;

/** @name Pool weight management fields
 * These fields manage the weights and probabilities of selected scenarios
 * @{ */

 /// Sum of the weights of scenarios in the discrete pool
 /** Holds the sum of weights for all selected scenarios in
  * \c scenarioIndexes.
  * Used as denominator when computing normalized probabilities. */
 double sumPoolWeights;
 
 /// Weights of scenarios in the selected pool
 /** Vector containing the weights of scenarios that have been
  * selected for the current pool (via init_random_pool or
  * init_representative_pool).
  * For random pool: renormalized original weights
  * For representative pool: aggregated weights from assignments
  * The size of this vector equals the number of selected scenarios.
  * This container is only populated after pool initialization. */
 std::vector< double > poolWeights;

/** @} ---------------------------------------------------------------------*/
/*-------------------- HELPER METHODS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name helper methods of the class
 * Miscellaneous functions
 * @{ */

 /// Empty the scenario pool
 /** Function to clear the internal pool of selected scenario indices.
  * Resets \c scenarioIndexes, \c sumPoolWeights, \c currentScenarioIndex,
  * \c poolSize, and \c is_initialized. */

 void empty_pool( void );

/*--------------------------------------------------------------------------*/
 /// Compute the transport cost matrix between scenarios
 /** Calculates the distance matrix between all pairs of scenarios using
  * the specified norm (ell parameter).
  * 
  * @param n_scenarios The total number of scenarios
  * @param scenario_size The dimension of each scenario
  * @param ell The power parameter for the distance calculation
  * @return The transport cost matrix
  */
 virtual CMatrix
 compute_transport_cost_matrix( ScenarioIndex n_scenarios ,
				ScenarioSize scenario_size ,
				float ell ) const;

/*--------------------------------------------------------------------------*/
 /// Compute distance between two scenarios
 /** Helper method to compute the ell-norm distance between two scenarios.
  * 
  * @param scenario1 First scenario as an Eigen vector
  * @param scenario2 Second scenario as an Eigen vector  
  * @param ell The power parameter for the distance calculation
  * @return The ell-power of the norm distance */

 virtual double compute_scenario_distance( const Eigen::VectorXd& scenario1 ,
                                           const Eigen::VectorXd& scenario2 ,
					   float ell ) const;

/*--------------------------------------------------------------------------*/
 /// Update pool weights after scenario selection (for random pool)
 /** Computes \c sumPoolWeights and populates \c poolWeights based on 
  * the selected scenarios in \c scenarioIndexes after \c init_random_pool().
  * 
  * This method:
  * 1. Calculates the sum of weights for all selected scenarios
  *    (\c sumPoolWeights)
  * 2. Populates the \c poolWeights container with normalized probabilities
  *    computed as: <tt>setWeights[scenario_i] / sumPoolWeights</tt>
  * 3. Falls back to uniform distribution if \c sumPoolWeights is zero.
  * 
  * @note This method is used by init_random_pool().
  *       For init_representative_pool(),
  *       use update_pool_weights_with_assignments() instead. */

 void update_pool_weights( void );
 
/*--------------------------------------------------------------------------*/
 /// Update pool weights with assignment information (for representative pool)
 /** Computes aggregated weights for representative scenarios based on the 
  * assignment of all scenarios to their closest representatives.
  * 
  * In scenario reduction, each selected scenario acts as a representative for
  * a cluster of scenarios. The weight of each representative is the sum of 
  * weights of all scenarios assigned to it (those for which it is the closest
  * selected scenario).
  * 
  * @param scenario_assignments Vector where element i contains the index of
  *                             the representative scenario that scenario i
  *                             is assigned to. Must have size nbScenarios.
  * 
  * @note This method is used by init_representative_pool() to properly
  *       aggregate weights based on the scenario reduction solution. */

 void update_pool_weights_with_assignments(
		const std::vector< ScenarioIndex > & scenario_assignments );
 
/*--------------------------------------------------------------------------*/
 /// Apply baseline selection method
 /** Selects the top scenarios based on their probability weights.
  * This is the fallback method when no BlockSolverConfig is provided.
  * Scenarios with higher weights are selected first.
  * 
  * @param target_size Number of scenarios to select */

 void apply_baseline_selection( ScenarioIndex target_size );

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };   // end( class DiscreteScenarioSet )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 } // end( namespace SMSpp_di_unipi_it )

#endif /* __DiscreteScenarioSet */

/*--------------------------------------------------------------------------*/
/*------------------- End file DiscreteScenarioSet.h -----------------------*/
/*--------------------------------------------------------------------------*/
