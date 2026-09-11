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

#include "DiscreteScenarioSet.h"

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
 v_pool.clear();
 f_self_view.reset();

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

 // no reduction to start with: the pool of a node is all of its children
 v_pool.resize( N );
 for( std::size_t n = 0 ; n < N ; ++n )
  reset_node_pool( static_cast< NodeIndex >( n ) );
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
/*------------------------- VIEW-BASED ACCESS ------------------------------*/
/*--------------------------------------------------------------------------*/

std::unique_ptr< MultiStageScenarioGenerator::View >
MultiStageDiscreteScenarioSet::root_view( void ) const
{
 if( f_nodes.empty() )
  throw( std::logic_error(
   "MultiStageDiscreteScenarioSet::root_view: empty tree" ) );
 return( std::make_unique< TreeView >( this , get_root() ) );
 }

/*--------------------------------------------------------------------------*/

MultiStageDiscreteScenarioSet::TreeView &
MultiStageDiscreteScenarioSet::self_view( void ) const
{
 if( f_nodes.empty() )
  throw( std::logic_error(
   "MultiStageDiscreteScenarioSet::self_view: empty tree" ) );
 if( ! f_self_view )
  f_self_view = std::make_unique< TreeView >( this , get_root() );
 return( * f_self_view );
 }

/*--------------------------------------------------------------------------*/
/*----------------- ScenarioGenerator FACE = THE ROOT VIEW -----------------*/
/*--------------------------------------------------------------------------*/
// seen as a plain (single-stage) ScenarioGenerator, this object is its own
// root view: everything is delegated to it

ScenarioGenerator::ScenarioIndex
MultiStageDiscreteScenarioSet::get_support_size( void )
{
 return( self_view().get_support_size() );
 }

/*--------------------------------------------------------------------------*/

ScenarioGenerator::ScenarioSize
MultiStageDiscreteScenarioSet::get_scenario_size( void ) const
{
 return( self_view().get_scenario_size() );
 }

/*--------------------------------------------------------------------------*/

ScenarioGenerator::Scenario
MultiStageDiscreteScenarioSet::get_current_scenario( void ) const
{
 return( self_view().get_current_scenario() );
 }

/*--------------------------------------------------------------------------*/

double MultiStageDiscreteScenarioSet::
                              get_current_scenario_probability( void ) const
{
 return( self_view().get_current_scenario_probability() );
 }

/*--------------------------------------------------------------------------*/

bool MultiStageDiscreteScenarioSet::next_scenario( void )
{
 return( self_view().next_scenario() );
 }

/*--------------------------------------------------------------------------*/

void MultiStageDiscreteScenarioSet::reset_pool( void )
{
 self_view().reset_pool();
 }

/*--------------------------------------------------------------------------*/

void MultiStageDiscreteScenarioSet::init_random_pool( ScenarioIndex )
{
 // prototype: the whole tree is the pool, nothing to (re)sample
 }

/*--------------------------------------------------------------------------*/

void MultiStageDiscreteScenarioSet::init_representative_pool(
						       ScenarioIndex size )
{
 self_view().init_representative_pool( size );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- POOL OF A NODE ---------------------------------*/
/*--------------------------------------------------------------------------*/

void MultiStageDiscreteScenarioSet::reset_node_pool( NodeIndex n ) const
{
 auto & pool = v_pool.at( n );
 pool.selected = f_nodes.at( n ).children;
 pool.weights.clear();
 pool.weights.reserve( pool.selected.size() );
 for( auto c : pool.selected )
  pool.weights.push_back( get_node_probability( c ) );
 }

/*--------------------------------------------------------------------------*/

void MultiStageDiscreteScenarioSet::reduce_node( NodeIndex n ,
						 ScenarioIndex size ) const
{
 auto & pool = v_pool.at( n );
 const auto & children = f_nodes.at( n ).children;

 // INFScenario, or asking for no fewer than there are, means "no
 // restriction": give the node back all of its children
 if( ( size == INFScenario ) || ( size >= children.size() ) ) {
  reset_node_pool( n );
  ++pool.generation;
  return;
  }

 if( size == 0 )
  throw( std::invalid_argument(
   "MultiStageDiscreteScenarioSet::reduce_node: the number of "
   "representatives must be positive" ) );

 if( children.empty() )
  throw( std::logic_error(
   "MultiStageDiscreteScenarioSet::reduce_node: node " +
   std::to_string( n ) + " is a leaf, it has no pool to reduce" ) );

 // the children of the node, with their conditional probabilities, are an
 // ordinary discrete scenario set: hand them to one, and let it reduce
 // them through a ScenarioReductionBlock and the configured Solver, so
 // that the single- and the multi-stage case share the algorithms
 DiscreteScenarioSet dss;
 std::vector< std::vector< double > > scenarios;
 std::vector< double > weights;
 scenarios.reserve( children.size() );
 weights.reserve( children.size() );
 for( auto c : children ) {
  const auto d = get_node_data( c );
  scenarios.emplace_back( d.begin() , d.end() );
  weights.push_back( get_node_probability( c ) );
  }
 dss.load_from_memory( scenarios , weights );

 if( f_solver_config )
  dss.set_solver_config( f_solver_config->clone() );
 if( f_stochastic_block )
  dss.set_Block( f_stochastic_block );

 dss.init_representative_pool( size );

 // read back which children have been selected, and with which
 // probabilities the discarded ones have been accumulated onto them
 const auto selected = dss.get_selected_scenarios();
 pool.selected.clear();
 pool.weights.clear();
 pool.selected.reserve( selected.size() );
 pool.weights.reserve( selected.size() );
 dss.reset_pool();
 for( std::size_t i = 0 ; i < selected.size() ; ++i ) {
  pool.selected.push_back( children.at( selected[ i ] ) );
  pool.weights.push_back( dss.get_current_scenario_probability() );
  if( i + 1 < selected.size() )
   dss.next_scenario();
  }

 ++pool.generation;
 }

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS TreeView ---------------------------------*/
/*--------------------------------------------------------------------------*/

const std::string &
MultiStageDiscreteScenarioSet::TreeView::private_name( void ) const
{
 static const std::string name( "MultiStageDiscreteScenarioSet::TreeView" );
 return( name );
 }

/*--------------------------------------------------------------------------*/

bool MultiStageDiscreteScenarioSet::TreeView::is_valid( void ) const
{
 // the pinned history is still the one this view was created on iff every
 // node of it has the generation it had then; walking it upwards from the
 // current node reads v_stamps backwards, root last
 auto it = v_stamps.rbegin();
 for( auto n = f_node ; n != InvalidNode ; n = f_parent->get_parent( n ) ) {
  if( ( it == v_stamps.rend() ) ||
      ( * it != f_parent->get_node_generation( n ) ) )
   return( false );
  ++it;
  }
 return( it == v_stamps.rend() );
 }

/*--------------------------------------------------------------------------*/

void MultiStageDiscreteScenarioSet::TreeView::check_valid(
					       const char * method ) const
{
 if( ! is_valid() )
  throw( std::logic_error(
   std::string( "MultiStageDiscreteScenarioSet::TreeView::" ) + method +
   ": the pinned history has been regenerated, this view is no longer "
   "valid" ) );
 }

/*--------------------------------------------------------------------------*/

ScenarioGenerator::ScenarioSize
MultiStageDiscreteScenarioSet::TreeView::get_scenario_size( void ) const
{
 check_valid( "get_scenario_size" );
 const auto & ch = f_parent->get_pool( f_node );
 if( ch.empty() )
  return( 0 );
 return( f_parent->get_node_data( ch.front() ).size() );
 }

/*--------------------------------------------------------------------------*/

ScenarioGenerator::ScenarioIndex
MultiStageDiscreteScenarioSet::TreeView::get_support_size( void )
{
 check_valid( "get_support_size" );
 return( static_cast< ScenarioIndex >(
                                f_parent->get_pool( f_node ).size() ) );
 }

/*--------------------------------------------------------------------------*/

ScenarioGenerator::Scenario
MultiStageDiscreteScenarioSet::TreeView::get_current_scenario( void ) const
{
 check_valid( "get_current_scenario" );
 const auto & ch = f_parent->get_pool( f_node );
 if( f_pos >= ch.size() )
  throw( std::out_of_range(
   "MultiStageDiscreteScenarioSet::TreeView::get_current_scenario: "
   "iteration past the last realization" ) );
 return( f_parent->get_node_data( ch[ f_pos ] ) );
 }

/*--------------------------------------------------------------------------*/

double MultiStageDiscreteScenarioSet::TreeView::
                              get_current_scenario_probability( void ) const
{
 check_valid( "get_current_scenario_probability" );
 const auto & ch = f_parent->get_pool( f_node );
 if( f_pos >= ch.size() )
  throw( std::out_of_range(
   "MultiStageDiscreteScenarioSet::TreeView::"
   "get_current_scenario_probability: iteration past the last "
   "realization" ) );
 return( f_parent->get_pool_weight( f_node , f_pos ) );
 }

/*--------------------------------------------------------------------------*/

bool MultiStageDiscreteScenarioSet::TreeView::next_scenario( void )
{
 check_valid( "next_scenario" );
 if( f_pos + 1 < f_parent->get_pool( f_node ).size() ) {
  ++f_pos;
  return( true );
  }
 return( false );
 }

/*--------------------------------------------------------------------------*/

bool MultiStageDiscreteScenarioSet::TreeView::descend( void )
{
 check_valid( "descend" );
 const auto & ch = f_parent->get_pool( f_node );
 if( f_pos >= ch.size() )
  return( false );                       // no realization selected
 const auto child = ch[ f_pos ];
 if( f_parent->get_pool( child ).empty() )
  return( false );                       // the selected one is a leaf
 f_node = child;                         // H <- ( H , x_t )
 f_pos = 0;
 v_stamps.push_back( f_parent->get_node_generation( child ) );
 return( true );
 }

/*--------------------------------------------------------------------------*/

bool MultiStageDiscreteScenarioSet::TreeView::climb( void )
{
 check_valid( "climb" );
 const auto p = f_parent->get_parent( f_node );
 if( p == InvalidNode )
  return( false );                       // already pinned at the root
 // re-select the realization this view had descended into, so that
 // climb() undoes descend()
 const auto & sibs = f_parent->get_pool( p );
 std::size_t i = 0;
 while( ( i < sibs.size() ) && ( sibs[ i ] != f_node ) )
  ++i;
 f_node = p;
 f_pos = i;
 v_stamps.pop_back();
 return( true );
 }

/*--------------------------------------------------------------------------*/

void MultiStageDiscreteScenarioSet::TreeView::init_representative_pool(
						       ScenarioIndex size )
{
 check_valid( "init_representative_pool" );
 f_parent->reduce_node( f_node , size );
 // the pool this view reads is the one it has just re-defined, hence the
 // view stays valid, unlike any other one pinned at this node or below it
 v_stamps.back() = f_parent->get_node_generation( f_node );
 f_pos = 0;
 }

/*--------------------------------------------------------------------------*/

std::unique_ptr< MultiStageScenarioGenerator::View >
MultiStageDiscreteScenarioSet::TreeView::clone( void ) const
{
 return( std::make_unique< TreeView >( * this ) );
 }

/*--------------------------------------------------------------------------*/

MultiStageScenarioGenerator::StageIndex
MultiStageDiscreteScenarioSet::TreeView::stage( void ) const
{
 check_valid( "stage" );
 return( f_parent->get_node_stage( f_node ) );
 }

/*--------------------------------------------------------------------------*/
/*--------------- End File MultiStageDiscreteScenarioSet.cpp ---------------*/
/*--------------------------------------------------------------------------*/
