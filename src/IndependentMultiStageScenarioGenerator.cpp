/*--------------------------------------------------------------------------*/
/*------ File IndependentMultiStageScenarioGenerator.cpp -------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the *concrete* class
 * IndependentMultiStageScenarioGenerator, an implementation of
 * MultiStageScenarioGenerator for the case where the random variables of
 * each stage are mutually independent.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "IndependentMultiStageScenarioGenerator.h"

#include <stdexcept>
#include <string>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------- FACTORY MANAGEMENT ----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( IndependentMultiStageScenarioGenerator );

/*--------------------------------------------------------------------------*/
/*-------------------- CONSTRUCTION / DESTRUCTION --------------------------*/
/*--------------------------------------------------------------------------*/

IndependentMultiStageScenarioGenerator::
               ~IndependentMultiStageScenarioGenerator() { clear_inners(); }

/*--------------------------------------------------------------------------*/

void IndependentMultiStageScenarioGenerator::clear_inners( void )
{
 for( auto * inner : inners )
  delete inner;
 inners.clear();
 current_stage = 0;
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- DESERIALIZE -----------------------------------*/
/*--------------------------------------------------------------------------*/

void IndependentMultiStageScenarioGenerator::deserialize(
                                            const netCDF::NcGroup & group )
{
 clear_inners();

 // Mandatory "NumberStages" dimension
 std::size_t T = 0;
 deserialize_dim( group , "NumberStages" , T , false );
 if( T == 0 )
  throw( std::invalid_argument(
   "IndependentMultiStageScenarioGenerator::deserialize: "
   "NumberStages must be positive" ) );

 inners.reserve( T );
 for( std::size_t t = 0 ; t < T ; ++t ) {
  const auto stage_name = "Stage_" + std::to_string( t );
  const auto stage_group = group.getGroup( stage_name );
  if( stage_group.isNull() ) {
   clear_inners();
   throw( std::invalid_argument(
    "IndependentMultiStageScenarioGenerator::deserialize: "
    "missing sub-group '" + stage_name + "'" ) );
   }

  auto * inner = ScenarioGenerator::new_ScenarioGenerator( stage_group );
  if( ! inner ) {
   clear_inners();
   throw( std::runtime_error(
    "IndependentMultiStageScenarioGenerator::deserialize: "
    "failed to construct inner :ScenarioGenerator for stage " +
    std::to_string( t ) ) );
   }

  if( dynamic_cast< MultiStageScenarioGenerator * >( inner ) ) {
   delete inner;
   clear_inners();
   throw( std::invalid_argument(
    "IndependentMultiStageScenarioGenerator::deserialize: "
    "inner generator for stage " + std::to_string( t ) +
    " is itself a :MultiStageScenarioGenerator, which is not "
    "allowed (nesting multi-stage generators is not supported)" ) );
   }

  inners.push_back( inner );
  }

 // After deserialization the outer is left in canonical walkable state:
 // current_stage = 0 and each inner is already in its own canonical
 // state (per its lazy-init contract).
 current_stage = 0;
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- SERIALIZE ------------------------------------*/
/*--------------------------------------------------------------------------*/

void IndependentMultiStageScenarioGenerator::serialize(
                                            netCDF::NcGroup & group ) const
{
 ScenarioGenerator::serialize( group );

 group.putAtt( "type" , "IndependentMultiStageScenarioGenerator" );
 group.addDim( "NumberStages" , inners.size() );

 for( std::size_t t = 0 ; t < inners.size() ; ++t ) {
  auto stage_group = group.addGroup( "Stage_" + std::to_string( t ) );
  inners[ t ]->serialize( stage_group );
  }
 }

/*--------------------------------------------------------------------------*/
/*----------------- METHODS INHERITED FROM ScenarioGenerator ---------------*/
/*--------------------------------------------------------------------------*/

void IndependentMultiStageScenarioGenerator::set_seed( unsigned long seed )
{
 // Derive distinct per-stage seeds from the master seed so that the
 // shuffles of different stages are not perfectly correlated. The
 // offset-by-t scheme is reproducible and cheap.
 for( std::size_t t = 0 ; t < inners.size() ; ++t )
  inners[ t ]->set_seed( seed + static_cast< unsigned long >( t ) );
 }

/*--------------------------------------------------------------------------*/

void IndependentMultiStageScenarioGenerator::set_Block( Block * block )
{
 for( auto * inner : inners )
  inner->set_Block( block );
 }

/*--------------------------------------------------------------------------*/

void IndependentMultiStageScenarioGenerator::set_config(
                                                   Configuration * config )
{
 // special case of the SimpleConfiguration< std::vector< Configuration * >
 if( auto vcfg =
     dynamic_cast< SimpleConfiguration< std::vector< Configuration * > >
                   * >( config ) ) {
  for( std::size_t i = 0 ; i < vcfg->f_value.size() ; ++i )
   if( i >= inners.size() )
    break;
   else
    inners[ i ]->set_config( vcfg->f_value[ i ] );

  return;
  }
 
 // Forward the same Configuration to every inner; each inner is free to
 // interpret or ignore it (default base impl is a no-op).
 for( auto * inner : inners )
  inner->set_config( config );
 }

/*--------------------------------------------------------------------------*/

ScenarioGenerator::ScenarioSize
IndependentMultiStageScenarioGenerator::get_scenario_size( void ) const
{
 if( inners.empty() )
  throw( std::logic_error(
   "IndependentMultiStageScenarioGenerator::get_scenario_size: "
   "no inner :ScenarioGenerator has been deserialized" ) );
 return( inners[ current_stage ]->get_scenario_size() );
 }

/*--------------------------------------------------------------------------*/

ScenarioGenerator::Scenario
IndependentMultiStageScenarioGenerator::get_current_scenario( void ) const
{
 if( inners.empty() )
  throw( std::logic_error(
   "IndependentMultiStageScenarioGenerator::get_current_scenario: "
   "no inner :ScenarioGenerator has been deserialized" ) );
 return( inners[ current_stage ]->get_current_scenario() );
 }

/*--------------------------------------------------------------------------*/

double IndependentMultiStageScenarioGenerator::
                              get_current_scenario_probability( void ) const
{
 if( inners.empty() )
  throw( std::logic_error(
   "IndependentMultiStageScenarioGenerator::"
   "get_current_scenario_probability: no inner has been deserialized" ) );
 return( inners[ current_stage ]->get_current_scenario_probability() );
 }

/*--------------------------------------------------------------------------*/

bool IndependentMultiStageScenarioGenerator::next_scenario( void )
{
 if( inners.empty() )
  return( false );
 return( inners[ current_stage ]->next_scenario() );
 }

/*--------------------------------------------------------------------------*/

void IndependentMultiStageScenarioGenerator::reset_pool( void )
{
 // Per :MultiStageScenarioGenerator semantics, reset_pool() acts on
 // the current stage only — it does not move the stage cursor and
 // does not touch the iteration of other stages. To also rewind the
 // stage cursor, call previous_stage( INFStage ).
 if( inners.empty() )
  throw( std::runtime_error(
   "IndependentMultiStageScenarioGenerator::reset_pool: "
   "no inner :ScenarioGenerator has been deserialized" ) );
 inners[ current_stage ]->reset_pool();
 }

/*--------------------------------------------------------------------------*/

bool IndependentMultiStageScenarioGenerator::is_pool_initialized( void ) const
{
 if( inners.empty() )
  return( false );
 for( const auto * inner : inners )
  if( ! inner->is_pool_initialized() )
   return( false );
 return( true );
 }

/*--------------------------------------------------------------------------*/

bool IndependentMultiStageScenarioGenerator::next_stage( void )
{
 if( current_stage + 1 >= inners.size() )
  return( false );
 ++current_stage;
 inners[ current_stage ]->reset_pool();
 return( true );
 }

/*--------------------------------------------------------------------------*/

void IndependentMultiStageScenarioGenerator::previous_stage(
                                                          StageIndex step )
{
 if( step == 0 )
  return;
 current_stage = ( step >= current_stage ) ? 0 : current_stage - step;
 inners[ current_stage ]->reset_pool();
 }

/*--------------------------------------------------------------------------*/
/*-------------------- POOL INITIALIZATION METHODS -------------------------*/
/*--------------------------------------------------------------------------*/
// init_*_pool( K ) acts on the inner of the *current stage* only, in
// keeping with the :MultiStageScenarioGenerator convention that all
// "current scenario" methods refer to the current stage. The caller
// initializes per-stage pools by looping with next_stage() and issuing
// a fresh init_*_pool( K_t ) at each stage.

void IndependentMultiStageScenarioGenerator::init_random_pool(
                                                         ScenarioIndex size )
{
 if( inners.empty() )
  throw( std::logic_error(
   "IndependentMultiStageScenarioGenerator::init_random_pool: "
   "no inner :ScenarioGenerator has been deserialized" ) );

 inners[ current_stage ]->init_random_pool( size );
 }

/*--------------------------------------------------------------------------*/

void IndependentMultiStageScenarioGenerator::init_representative_pool(
                                                         ScenarioIndex size )
{
 if( inners.empty() )
  throw( std::logic_error(
   "IndependentMultiStageScenarioGenerator::init_representative_pool: "
   "no inner :ScenarioGenerator has been deserialized" ) );

 inners[ current_stage ]->init_representative_pool( size );
 }

/*--------------------------------------------------------------------------*/
/*----- End File IndependentMultiStageScenarioGenerator.cpp ----------------*/
/*--------------------------------------------------------------------------*/
