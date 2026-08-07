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
 * walkable state: each inner already in its canonical full-universe
 * state (per the inner-specific lazy-init contract).
 *
 * ### Reading the scenarios: views
 *
 * As any :MultiStageScenarioGenerator, this class is read through
 * MultiStageScenarioGenerator::View, here implemented by StageView: since
 * the stages are mutually independent, a history H is irrelevant and a
 * position in the process boils down to the stage t it has reached, whose
 * pool is the inner of stage t. descend() / climb() move to t + 1 / t - 1,
 * and the generator itself, seen as a plain ScenarioGenerator, behaves
 * like the view pinned at stage 0.
 *
 * Note that, the cursor of each stage being held by the inner
 * :ScenarioGenerator of that stage, two StageView pinned at the *same*
 * stage do share it, and therefore cannot be used independently; views
 * pinned at distinct stages are instead completely independent, which is
 * all a stage-independent process needs.
 *
 * ### Pool semantics
 *
 * init_random_pool( K ) / init_representative_pool( K ) act on the inner
 * of the stage the view they are called on is pinned at, this is the
 * :MultiStageScenarioGenerator convention (see the comments on the base
 * class scalar forms in ScenarioGenerator.h). To (re-)init the pool of
 * more than one stage, walk them with View::descend() and call
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

 [[nodiscard]] ScenarioIndex get_support_size( void ) override;

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

 [[nodiscard]] bool is_stage_independent( void ) const override {
  return( true );
  }

/*--------------------------------------------------------------------------*/
 /// create a View pinned at stage 0 (the MSSG view entry point)
 /** Implements MultiStageScenarioGenerator::root_view(): returns a View (a
  * StageView) pinned at stage 0, whose pool is the inner of that stage.
  * The other stages are reached from it via View::descend(). */

 [[nodiscard]] std::unique_ptr< View > root_view( void ) const override;

/*--------------------------------------------------------------------------*/
 // init_random_pool() / init_representative_pool() act on the inner of
 // the stage of the view they are called on — see the comments on the
 // base class scalar forms in ScenarioGenerator.h. No per-stage or
 // vector overloads are exposed: the per-stage walk with View::descend()
 // is the canonical way to (re-)init more than one stage at a time.

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
/*--------------------------- PRIVATE TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
 /// read-only single-stage view of one stage of the process
 /** A StageView is the object returned by
  * IndependentMultiStageScenarioGenerator::root_view(). Since the stages
  * are mutually independent, the history H that a View pins is irrelevant
  * and only the stage t it has reached matters: the StageView is therefore
  * a thin forwarder to the inner :ScenarioGenerator of stage t, and
  * descend() / climb() just move to the next / previous stage. Because
  * that inner holds the cursor, two StageView pinned at the same stage
  * share it (see the class comments). */

 class StageView : public MultiStageScenarioGenerator::View
 {
  public:

   StageView( const IndependentMultiStageScenarioGenerator * parent ,
	      StageIndex stage )
    : f_parent( parent ) , f_stage( stage ) { }

   void deserialize( const netCDF::NcGroup & ) override {
    throw( std::logic_error(
     "IndependentMultiStageScenarioGenerator::StageView::deserialize: a "
     "view is not independently (de)serializable" ) );
    }

   void set_seed( unsigned long seed = 0 ) override {
    inner()->set_seed( seed );
    }

   [[nodiscard]] ScenarioIndex get_support_size( void ) override {
    return( inner()->get_support_size() );
    }

   [[nodiscard]] ScenarioSize get_scenario_size( void ) const override {
    return( inner()->get_scenario_size() );
    }

   [[nodiscard]] Scenario get_current_scenario( void ) const override {
    return( inner()->get_current_scenario() );
    }

   [[nodiscard]] double get_current_scenario_probability( void )
    const override { return( inner()->get_current_scenario_probability() ); }

   bool next_scenario( void ) override {
    return( inner()->next_scenario() );
    }

   void reset_pool( void ) override { inner()->reset_pool(); }

   [[nodiscard]] bool is_pool_initialized( void ) const override {
    return( inner()->is_pool_initialized() );
    }

   // the stages being independent, each is reduced on its own
   void init_random_pool( ScenarioIndex size = INFScenario ) override {
    inner()->init_random_pool( size );
    }

   void init_representative_pool( ScenarioIndex size = INFScenario )
    override { inner()->init_representative_pool( size ); }

   // history-pinned View interface: move through the stages
   bool descend( void ) override;
   bool climb( void ) override;
   [[nodiscard]] std::unique_ptr< View > clone( void ) const override;
   [[nodiscard]] StageIndex stage( void ) const override {
    return( f_stage );
    }

  private:

   /// the inner :ScenarioGenerator this view is pinned at
   [[nodiscard]] ScenarioGenerator * inner( void ) const;

   const IndependentMultiStageScenarioGenerator * f_parent;  ///< owner
   StageIndex f_stage;        ///< stage whose inner this view reads

   [[nodiscard]] const std::string & private_name( void ) const override;

 };   // end( class StageView )

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
