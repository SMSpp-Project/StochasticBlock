/*--------------------------------------------------------------------------*/
/*------------------ File ScenarioReductionBlock.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the ScenarioReductionBlock class.
 *
 * \author Minh Duc Pham \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Minh Duc Pham
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ScenarioReductionBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_1( ScenarioReductionBlock );

/*--------------------------------------------------------------------------*/
/*----------- METHODS OF ScenarioReductionBlock ----------------------------*/
/*--------------------------------------------------------------------------*/

void ScenarioReductionBlock::set_stochastic_block( Block * tssb )
{
 // Clear existing sub-Blocks first
 if( ! v_Block.empty() ) {
  v_Block.clear();
  }

 if( tssb ) {
  v_Block.push_back( tssb );
  // Do not call set_f_Block: we do not own this Block
  }
 }

/*--------------------------------------------------------------------------*/

void ScenarioReductionBlock::print( std::ostream & output ,
                                    char vlvl ) const
{
 output << std::endl << "ScenarioReductionBlock";

 if( f_scenario_generator )
  output << " with ScenarioGenerator ("
         << f_scenario_generator->classname() << ")";
 else
  output << " with no ScenarioGenerator";

 if( f_solution.is_set() )
  output << ", solution set with "
         << f_solution.selected_indices.size()
         << " representative scenarios";
 else
  output << ", no solution set yet";

 output << std::endl;
 }

/*--------------------------------------------------------------------------*/

void ScenarioReductionBlock::serialize( netCDF::NcGroup & group ) const
{
 Block::serialize( group );
 group.putAtt( "type" , "ScenarioReductionBlock" );

 // Solution is not serialized: it is produced at runtime by a Solver
 }

/*--------------------------------------------------------------------------*/

void ScenarioReductionBlock::deserialize( const netCDF::NcGroup & group )
{
 Block::deserialize( group );
 // Nothing to deserialize beyond the base class for now
 }

/*--------------------------------------------------------------------------*/
/*-------------- End File ScenarioReductionBlock.cpp -----------------------*/
/*--------------------------------------------------------------------------*/