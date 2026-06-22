/*--------------------------------------------------------------------------*/
/*----------------- File MultiStageDiscreteScenarioSet.cpp -----------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the *concrete* class MultiStageDiscreteScenarioSet, an
 * implementation of MultiStageScenarioGenerator for a general (history-
 * dependent) discrete multistage process given explicitly as a scenario
 * tree. Associated with the header file MultiStageDiscreteScenarioSet.h
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "MultiStageDiscreteScenarioSet.h"

#include <stdexcept>
#include <string>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------- FACTORY MANAGEMENT ----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( MultiStageDiscreteScenarioSet );

/*--------------------------------------------------------------------------*/
/*-------------------------- DESERIALIZE -----------------------------------*/
/*--------------------------------------------------------------------------*/

void MultiStageDiscreteScenarioSet::deserialize(
                                            const netCDF::NcGroup & group )
{
 f_nodes.clear();
 f_stage_size.clear();
 f_current_node = 0;

 // mandatory dimensions
 std::size_t T = 0 , N = 0 , D = 0;
 ::deserialize_dim( group , "NumberStages" , T , false );
 ::deserialize_dim( group , "NumberNodes" , N , false );
 ::deserialize_dim( group , "ScenarioDataSize" , D , false );

 if( ( T == 0 ) || ( N == 0 ) || ( D == 0 ) )
  throw( std::invalid_argument(
   "MultiStageDiscreteScenarioSet::deserialize: NumberStages, NumberNodes "
   "and ScenarioDataSize must all be positive" ) );

 f_number_stages = static_cast< StageIndex >( T );

 // per-stage realization sizes d_t
 {
  std::vector< unsigned int > ss( T );
  auto v = group.getVar( "StageScenarioSize" );
  if( v.isNull() )
   throw( std::invalid_argument(
    "MultiStageDiscreteScenarioSet::deserialize: missing variable "
    "'StageScenarioSize'" ) );
  v.getVar( ss.data() );
  f_stage_size.assign( ss.begin() , ss.end() );
  }

 // per-node fields
 std::vector< unsigned int > nstage( N ) , nparent( N );
 std::vector< double > nprob( N ) , ndata( N * D );

 auto vstage  = group.getVar( "NodeStage" );
 auto vparent = group.getVar( "NodeParent" );
 auto vprob   = group.getVar( "NodeProbability" );
 auto vdata   = group.getVar( "NodeData" );
 if( vstage.isNull() || vparent.isNull() || vprob.isNull() || vdata.isNull() )
  throw( std::invalid_argument(
   "MultiStageDiscreteScenarioSet::deserialize: one of the mandatory "
   "variables NodeStage / NodeParent / NodeProbability / NodeData is "
   "missing" ) );

 vstage.getVar( nstage.data() );
 vparent.getVar( nparent.data() );
 vprob.getVar( nprob.data() );
 vdata.getVar( ndata.data() );

 // build the nodes (no children yet)
 f_nodes.resize( N );
 for( std::size_t n = 0 ; n < N ; ++n ) {
  auto & nd = f_nodes[ n ];
  nd.stage = static_cast< StageIndex >( nstage[ n ] );
  nd.parent = ( nparent[ n ] >= N ) ? InvalidNode
                : static_cast< NodeIndex >( nparent[ n ] );
  nd.probability = nprob[ n ];
  if( nd.stage >= f_number_stages )
   throw( std::invalid_argument(
    "MultiStageDiscreteScenarioSet::deserialize: node " +
    std::to_string( n ) + " has stage out of range" ) );
  // copy only the d_{stage} meaningful entries of the padded row
  const auto d = f_stage_size[ nd.stage ];
  nd.data.assign( ndata.begin() + n * D , ndata.begin() + n * D + d );
  }

 // wire up children from the parent links (nodes are parent-before-child)
 for( std::size_t n = 0 ; n < N ; ++n ) {
  const auto p = f_nodes[ n ].parent;
  if( p == InvalidNode )
   continue;
  if( p >= n )
   throw( std::invalid_argument(
    "MultiStageDiscreteScenarioSet::deserialize: nodes are not listed "
    "parent-before-child (node " + std::to_string( n ) +
    " has parent " + std::to_string( p ) + ")" ) );
  f_nodes[ p ].children.push_back( static_cast< NodeIndex >( n ) );
  }

 f_current_node = get_root();
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- SERIALIZE ------------------------------------*/
/*--------------------------------------------------------------------------*/

void MultiStageDiscreteScenarioSet::serialize( netCDF::NcGroup & group ) const
{
 ScenarioGenerator::serialize( group );

 group.putAtt( "type" , "MultiStageDiscreteScenarioSet" );

 const std::size_t T = f_number_stages;
 const std::size_t N = f_nodes.size();
 std::size_t D = 0;
 for( auto d : f_stage_size )
  D = std::max( D , d );

 auto stageDim = group.addDim( "NumberStages" , T );
 auto nodeDim  = group.addDim( "NumberNodes" , N );
 auto dataDim  = group.addDim( "ScenarioDataSize" , D );

 // per-stage sizes
 std::vector< unsigned int > ss( f_stage_size.begin() , f_stage_size.end() );
 group.addVar( "StageScenarioSize" , netCDF::NcUint() , stageDim ).putVar(
                                                              ss.data() );

 // per-node fields
 std::vector< unsigned int > nstage( N ) , nparent( N );
 std::vector< double > nprob( N ) , ndata( N * D , 0.0 );
 for( std::size_t n = 0 ; n < N ; ++n ) {
  const auto & nd = f_nodes[ n ];
  nstage[ n ] = nd.stage;
  nparent[ n ] = ( nd.parent == InvalidNode )
                 ? static_cast< unsigned int >( N ) : nd.parent;
  nprob[ n ] = nd.probability;
  for( std::size_t k = 0 ; k < nd.data.size() ; ++k )
   ndata[ n * D + k ] = nd.data[ k ];
  }

 group.addVar( "NodeStage" , netCDF::NcUint() , nodeDim ).putVar(
                                                          nstage.data() );
 group.addVar( "NodeParent" , netCDF::NcUint() , nodeDim ).putVar(
                                                          nparent.data() );
 group.addVar( "NodeProbability" , netCDF::NcDouble() , nodeDim ).putVar(
                                                          nprob.data() );
 group.addVar( "NodeData" , netCDF::NcDouble() ,
               { nodeDim , dataDim } ).putVar( ndata.data() );
 }

/*--------------------------------------------------------------------------*/
/*------------------- VIEW-BASED (PARALLEL) ACCESS -------------------------*/
/*--------------------------------------------------------------------------*/

ScenarioGenerator * MultiStageDiscreteScenarioSet::make_view(
                                                       NodeIndex n ) const
{
 if( n >= f_nodes.size() )
  throw( std::out_of_range(
   "MultiStageDiscreteScenarioSet::make_view: node index out of range" ) );
 return( new NodeView( this , n ) );
 }

/*--------------------------------------------------------------------------*/

std::unique_ptr< MultiStageScenarioGenerator::View >
MultiStageDiscreteScenarioSet::root_view( void ) const
{
 return( std::make_unique< NodeView >( this , get_root() ) );
 }

/*--------------------------------------------------------------------------*/
/*----------------- LEGACY SINGLE-CURSOR ACCESS ----------------------------*/
/*--------------------------------------------------------------------------*/

ScenarioGenerator::ScenarioSize
MultiStageDiscreteScenarioSet::get_scenario_size( void ) const
{
 if( f_nodes.empty() )
  throw( std::logic_error(
   "MultiStageDiscreteScenarioSet::get_scenario_size: empty tree" ) );
 return( f_nodes[ f_current_node ].data.size() );
 }

/*--------------------------------------------------------------------------*/

ScenarioGenerator::Scenario
MultiStageDiscreteScenarioSet::get_current_scenario( void ) const
{
 if( f_nodes.empty() )
  throw( std::logic_error(
   "MultiStageDiscreteScenarioSet::get_current_scenario: empty tree" ) );
 const auto & d = f_nodes[ f_current_node ].data;
 return( Scenario( d.data() , d.size() ) );
 }

/*--------------------------------------------------------------------------*/

double MultiStageDiscreteScenarioSet::
                              get_current_scenario_probability( void ) const
{
 if( f_nodes.empty() )
  throw( std::logic_error(
   "MultiStageDiscreteScenarioSet::get_current_scenario_probability: "
   "empty tree" ) );
 return( f_nodes[ f_current_node ].probability );
 }

/*--------------------------------------------------------------------------*/

bool MultiStageDiscreteScenarioSet::next_scenario( void )
{
 // another realization of X_t with the same history = next sibling
 const auto p = f_nodes[ f_current_node ].parent;
 if( p == InvalidNode )
  return( false );   // the root is unique
 const auto & sibs = f_nodes[ p ].children;
 for( std::size_t i = 0 ; i < sibs.size() ; ++i )
  if( sibs[ i ] == f_current_node ) {
   if( i + 1 < sibs.size() ) {
    f_current_node = sibs[ i + 1 ];
    return( true );
    }
   return( false );
   }
 return( false );
 }

/*--------------------------------------------------------------------------*/

void MultiStageDiscreteScenarioSet::reset_pool( void )
{
 // rewind to the first sibling at the current node's history
 const auto p = f_nodes.at( f_current_node ).parent;
 if( p != InvalidNode )
  f_current_node = f_nodes[ p ].children.front();
 }

/*--------------------------------------------------------------------------*/

bool MultiStageDiscreteScenarioSet::next_stage( void )
{
 // descend to the first child = first realization of X_{t+1} | history
 const auto & ch = f_nodes[ f_current_node ].children;
 if( ch.empty() )
  return( false );
 f_current_node = ch.front();
 return( true );
 }

/*--------------------------------------------------------------------------*/

void MultiStageDiscreteScenarioSet::previous_stage( StageIndex step )
{
 for( StageIndex s = 0 ; s < step ; ++s ) {
  const auto p = f_nodes[ f_current_node ].parent;
  if( p == InvalidNode )
   break;
  f_current_node = p;
  }
 }

/*--------------------------------------------------------------------------*/

void MultiStageDiscreteScenarioSet::init_random_pool( ScenarioIndex )
{
 // prototype: the whole tree is the pool, nothing to (re)sample
 }

/*--------------------------------------------------------------------------*/

void MultiStageDiscreteScenarioSet::init_representative_pool( ScenarioIndex )
{
 // prototype: identity reduction. This is where GLOBAL scenario-tree
 // reduction would rebuild a smaller equivalent tree; see the class
 // comments in the header.
 }

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS NodeView ---------------------------------*/
/*--------------------------------------------------------------------------*/

const std::string &
MultiStageDiscreteScenarioSet::NodeView::private_name( void ) const
{
 static const std::string name( "MultiStageDiscreteScenarioSet::NodeView" );
 return( name );
 }

/*--------------------------------------------------------------------------*/

ScenarioGenerator::ScenarioSize
MultiStageDiscreteScenarioSet::NodeView::get_scenario_size( void ) const
{
 const auto & ch = f_parent->f_nodes.at( f_node ).children;
 if( ch.empty() )
  return( 0 );
 return( f_parent->f_nodes[ ch.front() ].data.size() );
 }

/*--------------------------------------------------------------------------*/

ScenarioGenerator::ScenarioIndex
MultiStageDiscreteScenarioSet::NodeView::get_support_size( void )
{
 return( static_cast< ScenarioIndex >(
          f_parent->f_nodes.at( f_node ).children.size() ) );
 }

/*--------------------------------------------------------------------------*/

ScenarioGenerator::Scenario
MultiStageDiscreteScenarioSet::NodeView::get_current_scenario( void ) const
{
 const auto & ch = f_parent->f_nodes.at( f_node ).children;
 if( f_pos >= ch.size() )
  throw( std::out_of_range(
   "MultiStageDiscreteScenarioSet::NodeView::get_current_scenario: "
   "iteration past the last child" ) );
 const auto & d = f_parent->f_nodes[ ch[ f_pos ] ].data;
 return( Scenario( d.data() , d.size() ) );
 }

/*--------------------------------------------------------------------------*/

double MultiStageDiscreteScenarioSet::NodeView::
                              get_current_scenario_probability( void ) const
{
 const auto & ch = f_parent->f_nodes.at( f_node ).children;
 if( f_pos >= ch.size() )
  throw( std::out_of_range(
   "MultiStageDiscreteScenarioSet::NodeView::"
   "get_current_scenario_probability: iteration past the last child" ) );
 return( f_parent->f_nodes[ ch[ f_pos ] ].probability );
 }

/*--------------------------------------------------------------------------*/

bool MultiStageDiscreteScenarioSet::NodeView::next_scenario( void )
{
 const auto & ch = f_parent->f_nodes.at( f_node ).children;
 if( f_pos + 1 < ch.size() ) {
  ++f_pos;
  return( true );
  }
 return( false );
 }

/*--------------------------------------------------------------------------*/

std::unique_ptr< MultiStageScenarioGenerator::View >
MultiStageDiscreteScenarioSet::NodeView::descend( void ) const
{
 const auto & ch = f_parent->f_nodes.at( f_node ).children;
 if( f_pos >= ch.size() )
  return( nullptr );                       // no current child
 const auto child = ch[ f_pos ];
 if( f_parent->f_nodes[ child ].children.empty() )
  return( nullptr );                       // current child is a leaf
 return( std::make_unique< NodeView >( f_parent , child ) );
 }

/*--------------------------------------------------------------------------*/

MultiStageScenarioGenerator::StageIndex
MultiStageDiscreteScenarioSet::NodeView::stage( void ) const
{
 return( f_parent->f_nodes.at( f_node ).stage );
 }

/*--------------------------------------------------------------------------*/
/*--------------- End File MultiStageDiscreteScenarioSet.cpp ---------------*/
/*--------------------------------------------------------------------------*/
