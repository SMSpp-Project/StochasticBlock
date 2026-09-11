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

#ifndef __DiscreteScenarioSet
#define __DiscreteScenarioSet

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"
#include "BlockSolverConfig.h"
#include "ScenarioGenerator.h"
#include "ScenarioReductionBlock.h"

#include <Eigen/Dense>
#include <random>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS DiscreteScenarioSet --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
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
 * DiscreteScenarioSet provides two primary methods for selecting scenario
 * subsets:
 *
 * - **Random Selection** (init_random_pool()): randomly samples poolsize
 *   scenarios from the full set, always available without additional
 *   configuration.
 *
 * - **Scenario Reduction** (init_representative_pool()): selects a
 *   representative subset. This requires an appropriate
 *   Configuration, that can be passed either via set_config() or loaded
 *   during netCDF deserialization.
 *
 * ### Scenario Reduction
 *
 * When init_representative_pool(K) is called with a finite K, the class:
 *  1. Creates a ScenarioReductionBlock holding a pointer to itself.
 *  2. If a TwoStageStochasticBlock has been registered via set_Block(),
 *     adds it as the first sub-Block of the ScenarioReductionBlock so that
 *     CSSCScenarioReductionSolver can access the scenario sub-problems.
 *  3. Attaches the registered Solver (set via set_solver_config()) and
 *     calls compute().
 *  4. Reads the ScenarioReductionBlockSolution and updates the pool.
*/

class DiscreteScenarioSet : public ScenarioGenerator {

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/

 using DiscreteScenarioPool = boost::multi_array< double , 2 >;
 using CMatrix              = boost::multi_array< double , 2 >;

/*--------------------------------------------------------------------------*/
/*---------- CONSTRUCTING AND DESTRUCTING DiscreteScenarioSet --------------*/
/*--------------------------------------------------------------------------*/

 DiscreteScenarioSet();

 ~DiscreteScenarioSet() override;

 void deserialize( const netCDF::NcGroup & group ) override;

 void serialize( netCDF::NcGroup & group ) const override;

 /** Populate from in-memory data. */
 void load_from_memory(
  const std::vector< std::vector< double > > & scenarios ,
  const std::vector< double > & weights = {} );

/*--------------------------------------------------------------------------*/
/*--------------------- ScenarioGenerator INTERFACE -----------------------*/
/*--------------------------------------------------------------------------*/

 void set_seed( unsigned long seed ) override { rng.seed( seed ); }

 [[nodiscard]] inline Scenario get_current_scenario() const override {
  const auto index = scenarioIndexes[ currentScenarioIndex ];
  return Scenario( &scenarioSet[ index ][ 0 ] , scenarioSize );
  }

 [[nodiscard]] double get_current_scenario_probability() const override {
  return poolWeights[ currentScenarioIndex ];
  }

 bool next_scenario() override;

 [[nodiscard]] ScenarioSize get_scenario_size() const override {
  return scenarioSize;
  }

 [[nodiscard]] ScenarioIndex get_support_size() override {
  return nbScenarios;
  }

 /** Set the partner Block. If block is a TwoStageStochasticBlock, stores
  * it so that init_representative_pool() can pass it to the
  * ScenarioReductionBlock. The default ScenarioGenerator implementation
  * does nothing; we override here to capture the pointer. */
 void set_Block( Block * block ) override {
  f_stochastic_block = block;
  }

/*--------------------------------------------------------------------------*/
/*----------------------- POOL MANAGEMENT ---------------------------------*/
/*--------------------------------------------------------------------------*/

 void init_random_pool( ScenarioIndex pool_size = INFScenario ) override;

 /** Select K representative scenarios via a ScenarioReductionBlock.
  * Requires a Solver to have been registered via set_solver_config().
  * Falls back to baseline (highest-weight) selection if no Solver is set.
  * @param target_pool_size  K, number of representatives to select. */
 void init_representative_pool(
  ScenarioIndex target_pool_size = INFScenario ) override;

/*--------------------------------------------------------------------------*/
/*----------------------- CONFIGURATION -----------------------------------*/
/*--------------------------------------------------------------------------*/

 /** Set the BlockSolverConfig used to attach a Solver to the
  * ScenarioReductionBlock inside init_representative_pool(). Ownership is
  * transferred. Resets f_solution_ready. */
 void set_solver_config( BlockSolverConfig * solver_config ) {
  f_solver_config.reset( solver_config );
  f_solution_ready = false;
  }

 /** Convenience overload also setting the pool size. Resets
  * f_solution_ready. */
 void set_solver_config( BlockSolverConfig * solver_config ,
                         ScenarioIndex pool_size ) {
  f_solver_config.reset( solver_config );
  poolSize = pool_size;
  f_solution_ready = false;
  }

 /** Mark that a valid solution is available. Called by
  * init_representative_pool() after the Solver writes its solution into
  * the ScenarioReductionBlock. */
 void set_solution_ready() { f_solution_ready = true; }

 /** Set configuration from a generic Configuration object.
  * Supported: SimpleConfiguration<int> (baseline), or
  * SimpleConfiguration<pair<int,Configuration*>> (with solver). */
 void set_config( Configuration * config ) override;

/*--------------------------------------------------------------------------*/
/*------------------------------ GETTERS ----------------------------------*/
/*--------------------------------------------------------------------------*/

 [[nodiscard]] const ScenarioIndex & get_nbScenarios() const {
  return nbScenarios;
  }

 [[nodiscard]] const ScenarioSize & get_scenarioSize() const {
  return scenarioSize;
  }

 [[nodiscard]] bool is_pool_initialized() const override {
  return is_initialized;
  }

 void reset_pool() override {
  if( ! is_initialized )
   throw std::runtime_error( "reset_pool: pool not initialized" );
  currentScenarioIndex = 0;
  }

 [[nodiscard]] double get_scenario_value(
  ScenarioIndex scenario_idx , ScenarioSize component_idx ) const {
  if( scenario_idx >= nbScenarios || component_idx >= scenarioSize )
   throw std::out_of_range( "get_scenario_value: index out of range" );
  return scenarioSet[ scenario_idx ][ component_idx ];
  }

 [[nodiscard]] std::span< const ScenarioIndex > get_selected_scenarios()
  const;

 [[nodiscard]] ScenarioIndex get_poolSize() const { return poolSize; }

 [[nodiscard]] float get_ell() const { return ell; }

 void set_ell( float ell_value ) {
  if( ell_value <= 0 )
   throw std::invalid_argument( "set_ell: ell must be positive" );
  ell = ell_value;
  }

 [[nodiscard]] inline std::span< const double > get_set_weights() const {
  return setWeights;
  }

 [[nodiscard]] inline std::span< const double > get_pool_weights() const {
  if( ! is_initialized )
   throw std::runtime_error( "get_pool_weights: pool not initialized" );
  return poolWeights;
  }

 /** Returns the number of representative scenarios that
  * init_representative_pool() will select. Used by
  * ScenarioReductionSolver to determine K. */
 [[nodiscard]] ScenarioIndex get_NReducedScenarios() const {
  return poolSize;
  }

 /** Returns the full weight vector (all N scenarios). Used by
  * ScenarioReductionSolver to build the internal CFLB. */
 [[nodiscard]] const std::vector< double > & get_Weights() const {
  return setWeights;
  }

 /** Returns a single scenario as a vector. Used by
  * ScenarioReductionSolver to compute pairwise distances.
  * @param idx  Scenario index in [0, nbScenarios). */
 [[nodiscard]] std::vector< double > get_scenario(
  ScenarioIndex idx ) const {
  if( idx >= nbScenarios )
   throw std::out_of_range( "get_scenario: index out of range" );
  return std::vector< double >(
   &scenarioSet[ idx ][ 0 ] ,
   &scenarioSet[ idx ][ 0 ] + scenarioSize );
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 DiscreteScenarioPool scenarioSet;
 std::vector< ScenarioIndex > scenarioIndexes;
 std::vector< ScenarioIndex > universeIndexes;

 ScenarioIndex currentScenarioIndex{ 0 };
 ScenarioIndex nbScenarios  = 0;
 ScenarioSize  scenarioSize = 0;
 ScenarioIndex poolSize     = 0;

 std::mt19937 rng;
 bool is_initialized{ false };

 float ell = 2.0f;  ///< power for the Wasserstein distance

 /** Pointer to the TwoStageStochasticBlock set by set_Block(). Not owned. */
 Block * f_stochastic_block = nullptr;

 /** Solver configuration for scenario reduction. Owned. */
 std::unique_ptr< BlockSolverConfig > f_solver_config;

 /** True when a valid solution has been written into the
  * ScenarioReductionBlock. init_representative_pool() uses the existing
  * solution instead of invoking the Solver again. Reset to false whenever
  * anything changes that would invalidate the solution (new data, new K,
  * new solver config). */
 bool f_solution_ready = false;

 std::vector< double > setWeights;
 double sumPoolWeights = 0.0;
 std::vector< double > poolWeights;

/*--------------------------------------------------------------------------*/
/*------------------------ PRIVATE METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

 void empty_pool();

 virtual CMatrix compute_transport_cost_matrix(
  ScenarioIndex n_scenarios , ScenarioSize scenario_size , float ell ) const;

 virtual double compute_scenario_distance(
  const Eigen::VectorXd & s1 , const Eigen::VectorXd & s2 , float ell ) const;

 void update_pool_weights();

 void update_pool_weights_with_assignments(
  const std::vector< ScenarioIndex > & assignments );

 void apply_baseline_selection( ScenarioIndex target_size );

 SMSpp_insert_in_factory_h;

 };  // end( class DiscreteScenarioSet )

}  // end( namespace SMSpp_di_unipi_it )

#endif  /* __DiscreteScenarioSet */

/*--------------------------------------------------------------------------*/
/*------------------- End file DiscreteScenarioSet.h -----------------------*/
/*--------------------------------------------------------------------------*/