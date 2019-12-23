/*--------------------------------------------------------------------------*/
/*------------------------ File StochasticBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the StochasticBlock class.
 *
 * \version 0.10
 *
 * \date 09 - 12 - 2019
 *
 * \author Rafael Durbano Lobato \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy; by Rafael Durbano Lobato
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

void StochasticBlock::deserialize( netCDF::NcGroup & group ) {
 auto inner_block_group = group.getGroup( "Block" );
 if( inner_block_group.isNull() )
  throw std::invalid_argument( "StochasticBlock::deserialize: the 'Block' "
                               "group must be present in the given "
                               "netCDF::NcGroup." );

 auto inner_block = new_Block( inner_block_group, this );
 if( ! inner_block )
  throw std::logic_error( "StochasticBlock::deserialize: the 'Block'"
                          "group is present but its description is "
                          "incomplete." );

 set_inner_block( inner_block );

 Index num_data_mappings;
 ::deserialize( group , "NumDataMappings" , & num_data_mappings , false );

 data_mappings.clear();
 data_mappings.reserve( num_data_mappings );

 for( Index i = 0 ; i < num_data_mappings ; ++i ) {
  auto data_mapping_group = group.getGroup( "DataMapping_" +
                                            std::to_string( i ) );

  if( data_mapping_group.isNull() )
   throw std::logic_error( "StochasticBlock::deserialize: 'DataMapping_" +
                           std::to_string( i ) + "' sub-group must be "
                           "present." );

  std::string template_parameter_types;
  if( ! ::deserialize( data_mapping_group , "TemplateParameterTypes" ,
                       & template_parameter_types , true ) )
   data_mappings.emplace_back( new SimpleDataMapping<> );
  else
   data_mappings.emplace_back( SimpleDataMappingFactory::new_SimpleDataMapping
                               ( template_parameter_types ) );
  data_mappings.back()->deserialize( data_mapping_group , inner_block );
 }
}

/*--------------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/

void StochasticBlock::add_Modification( sp_Mod mod ,
                                        Observer::ChnlName chnl ) {
 // TODO
 if( anyone_there() )
  add_Modification( std::make_shared<NBModification>( this ) );
}

/*--------------------------------------------------------------------------*/
/*------- METHODS FOR LOADING, PRINTING & SAVING THE StochasticBlock -------*/
/*--------------------------------------------------------------------------*/

void StochasticBlock::serialize( netCDF::NcGroup & group ) const {

 group.putAtt( "type" , "StochasticBlock" );

 auto inner_block = get_inner_block();

 if( inner_block ) {
  auto inner_block_group = group.addGroup( "Block" );
  v_Block[ 0 ]->serialize( inner_block_group );
 }

 ::serialize( group , "NumDataMappings" , netCDF::NcUint64() ,
              data_mappings.size() );

 for( Index i = 0 ; i < data_mappings.size() ; ++i ) {
  auto data_mapping_group = group.addGroup( "DataMapping_" +
                                            std::to_string( i ) );
  data_mappings[ i ]->serialize( data_mapping_group , inner_block );
 }
}

/*--------------------------------------------------------------------------*/

void StochasticBlock::print( std::ostream &output ) const {
 output << std::endl << "StochasticBlock with ";

 if( v_Block.empty() )
  output << "no inner Block";
 else
  output << "the inner Block " << v_Block[ 0 ] << std::endl;
}

/*--------------------------------------------------------------------------*/
/*-------------------- End File StochasticBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
