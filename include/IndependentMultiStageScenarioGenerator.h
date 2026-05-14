/*--------------------------------------------------------------------------*/
/*------- File IndependentMultiStageScenarioGenerator.h --------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class
 * IndependentMultiStageScenarioGenerator, an implementation of
 * MultiStageScenarioGenerator for the case where the random variables of
 * each stage are mutually independent. The class is just a vector of
 * (single-stage) ScenarioGenerator *, one for each stage, and all the
 * multi-stage operations dispatch to the corresponding inner generator.
 *
 * No assumption is made on the concrete type of the inner generators: any
 * :ScenarioGenerator (registered in the factory) can be used as the inner
 * of a given stage. The only restriction is that an inner cannot itself
 * be a :MultiStageScenarioGenerator.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Claude Opus 4.7 \n
 *         Antrophic \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __IndependentMultiStageScenarioGenerator
 #define __IndependentMultiStageScenarioGenerator
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ScenarioGenerator.h"

#include <vector>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*-------- CLASS IndependentMultiStageScenarioGenerator --------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// concrete MultiStageScenarioGenerator for stage-independent processes
/** IndependentMultiStageScenarioGenerator implements the case where the
 * random variables ( X_0 , ... , X_{T-1} ) of the multistage process are
 * mutually independent: each X_t can be drawn without any reference to
 * the history H_t. The class is just a vector of T single-stage
 * :ScenarioGenerator pointers, one for each stage, and all multi-stage
 * operations dispatch to the inner of the current stage.
 *
 * ### Inner generators
 *
 * Any :ScenarioGenerator can be used as an inner, with the single
 * exception that an inner cannot itself be a :MultiStageScenarioGenerator
 * (the contract of this class is that each stage exposes a single
 * scalar random variable, so nesting multi-stage generators is not
 * meaningful). Inners are owned by this object: they are constructed in
 * deserialize() via the ScenarioGenerator factory and destroyed in the
 * destructor.
 *
 * ### netCDF format
 *
 * A group containing an IndependentMultiStageScenarioGenerator has a
 * mandatory string attribute "type" = "IndependentMultiStageScenario
 * Generator" plus:
 *
 *  - dimension "NumberStages" (size T > 0);
 *
 *  - one subgroup "Stage_0", ..., "Stage_{T-1}", each containing a
 *    :ScenarioGenerator (in factory format, i.e., with its own "type"
 *    attribute, or as an indirect group with "filename" attribute).
 *
 * After deserialize() returns, the generator is left in canonical
 * walkable state: current_stage = 0, each inner already in its canonical
 * full-universe state (per the inner-specific lazy-init contract).
 *
 * ### Pool semantics
 *
 * init_random_pool( K ) / init_representative_pool( K ) act on the
 * inner of the *current stage* only — this is the
 * :MultiStageScenarioGenerator convention (see the comments on the
 * base class scalar forms in ScenarioGenerator.h). To (re-)init the
 * pool of more than one stage, walk them with next_stage() and call
 * init_*_pool() once per stage with the desired per-stage size. */

class IndependentMultiStageScenarioGenerator
      : public MultiStageScenarioGenerator
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--- CONSTRUCTING AND DESTRUCTING IndependentMultiStageScenarioGenerator --*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing IndependentMultiStageScenarioGenerator
 *  @{ */

 IndependentMultiStageScenarioGenerator( void )
  : MultiStageScenarioGenerator() { }

/*--------------------------------------------------------------------------*/

 virtual ~IndependentMultiStageScenarioGenerator();

/*--------------------------------------------------------------------------*/

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*------- METHODS INHERITED FROM (Multi-Stage)ScenarioGenerator ------------*/
/*--------------------------------------------------------------------------*/
/** @name Functions inherited from (Multi-Stage)ScenarioGenerator
 *  @{ */

 void set_seed( unsigned long seed = 0 ) override;

/*--------------------------------------------------------------------------*/

 void set_Block( Block * block ) override;

/*--------------------------------------------------------------------------*/
 /// set the Configuration for the IndependentMultiStageScenarioGenerator
 /** The Configuration of the IndependentMultiStageScenarioGenerator can be
  * of two different types:
  *
  * - a SimpleConfiguration< std::vector< Configuration * >, in which case
  *   its i-th element is used as the Configuration for the ScenarioGenerator
  *   of the i-th stage (if the vector is shorter than the number of stages
  *   the uncovered stages are left un-Configure-d, if it is longer the extra
  *   entries are ignored);
  *
  * - any other type of Configuration, in which case it is passed verbatim
  *   to all the ScenarioGenerator of every stage (be careful, then, if it
  *   is / contains a BlockConfig as it gets "consumed" when apply()-ed). */

 void set_config( Configuration * config ) override;

/*--------------------------------------------------------------------------*/

 [[nodiscard]] ScenarioSize get_scenario_size( void ) const override;

/*--------------------------------------------------------------------------*/

 [[nodiscard]] Scenario get_current_scenario( void ) const override;

/*--------------------------------------------------------------------------*/

 [[nodiscard]] double get_current_scenario_probability( void ) const override;

/*--------------------------------------------------------------------------*/

 bool next_scenario( void ) override;

/*--------------------------------------------------------------------------*/

 void reset_pool( void ) override;

/*--------------------------------------------------------------------------*/

 [[nodiscard]] bool is_pool_initialized( void ) const override;

/*--------------------------------------------------------------------------*/

 [[nodiscard]] StageIndex get_stage_number( void ) override {
  return( static_cast< StageIndex >( inners.size() ) );
  }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] StageIndex get_current_stage( void ) override {
  return( current_stage );
  }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] bool next_stage( void ) override;

/*--------------------------------------------------------------------------*/

 void previous_stage( StageIndex step = 1 ) override;

/*--------------------------------------------------------------------------*/

 [[nodiscard]] bool is_stage_independent( void ) const override {
  return( true );
  }

/*--------------------------------------------------------------------------*/
 // init_random_pool() / init_representative_pool() act on the inner of
 // the current stage only — see the comments on the base class scalar
 // forms in ScenarioGenerator.h. No per-stage or vector overloads are
 // exposed: the per-stage loop with next_stage() is the canonical way
 // to (re-)init more than one stage at a time.

 void init_random_pool( ScenarioIndex size = INFScenario ) override;

 void init_representative_pool( ScenarioIndex size = INFScenario ) override;

/** @} ---------------------------------------------------------------------*/
/*----------------------- GETTERS / DIRECT ACCESS --------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Getters / direct access to inners
 *  @{ */

 /// returns the inner :ScenarioGenerator of stage \p t (read-only)
 [[nodiscard]] const ScenarioGenerator * get_inner( StageIndex t ) const {
  return( inners.at( t ) );
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

 /// release all the inner :ScenarioGenerator objects
 void clear_inners( void );

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 /// inner per-stage :ScenarioGenerator pointers, owned by this object
 std::vector< ScenarioGenerator * > inners;

 /// current stage in 0, ..., inners.size() - 1
 StageIndex current_stage = 0;

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

};   // end( class IndependentMultiStageScenarioGenerator )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

#endif  /* IndependentMultiStageScenarioGenerator.h included */

/*--------------------------------------------------------------------------*/
/*------ End File IndependentMultiStageScenarioGenerator.h -----------------*/
/*--------------------------------------------------------------------------*/
