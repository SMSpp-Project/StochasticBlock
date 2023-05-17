/*--------------------------------------------------------------------------*/
/*------------------------ File StochasticBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the StochasticBlock class.
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Rafael Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "StochasticBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_1( StochasticBlock );

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS of StochasticBlock ----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------- CONSTRUCTING AND DESTRUCTING StochasticBlock ---------------*/
/*--------------------------------------------------------------------------*/

void StochasticBlock::deserialize( const netCDF::NcGroup & group ) {
 auto inner_block_group = group.getGroup( "Block" );
 if( ! inner_block_group.isNull() ) {

  auto inner_block = new_Block( inner_block_group, this );
  if( ! inner_block )
   throw std::logic_error( "StochasticBlock::deserialize: the 'Block' "
                           "group is present but its description is "
                           "incomplete." );
  set_inner_block( inner_block );
 }

 data_mappings.clear();
 Index num_data_mappings;
 if( ::SMSpp_di_unipi_it::deserialize_dim( group , "NumberDataMappings" ,
                                           num_data_mappings , true ) &&
     num_data_mappings > 0 ) {
  data_mappings.reserve( num_data_mappings );
  assert( v_Block.size() == 1 && v_Block.front() );
  SimpleDataMappingBase::deserialize( group , data_mappings , v_Block.front() );
 }

 Block::deserialize( group );
}

/*--------------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/

void StochasticBlock::add_Modification( sp_Mod mod ,
                                        Observer::ChnlName chnl ) {
 // TODO
 if( anyone_there() )
  Block::add_Modification( std::make_shared< NBModification >( this ) );
}

/*--------------------------------------------------------------------------*/
/*----------- METHODS FOR PRINTING & SAVING THE StochasticBlock ------------*/
/*--------------------------------------------------------------------------*/

void StochasticBlock::print( std::ostream & output , char vlvl ) const
{
 output << std::endl << "StochasticBlock with ";

 if( v_Block.empty() )
  output << "no inner Block";
 else
  output << "the inner Block " << v_Block.front() << std::endl;
 }

/*--------------------------------------------------------------------------*/

void StochasticBlock::serialize( netCDF::NcGroup & group ) const
{
 Block::serialize( group );

 group.putAtt( "type" , "StochasticBlock" );

 Block * inner_block = nullptr;
 if( ! v_Block.empty() ) {
  assert( v_Block.size() == 1 );
  inner_block = v_Block.front();
  }

 if( inner_block ) {
  auto inner_block_group = group.addGroup( "Block" );
  inner_block->serialize( inner_block_group );
  }

 SimpleDataMappingBase::serialize( group , data_mappings , inner_block );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- End File StochasticBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
