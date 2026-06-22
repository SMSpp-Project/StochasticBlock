/*--------------------------------------------------------------------------*/
/*-------------------------- File test_mssds.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Stand-alone test / demo for MultiStageDiscreteScenarioSet, the discrete
 * (history-dependent) scenario-tree MultiStageScenarioGenerator.
 *
 * It builds a small 3-stage scenario tree
 *
 *                       root (stage 0)
 *                      /              \
 *                 s0 (0.6)          s1 (0.4)          <- stage 1
 *               /   |   \          /   |   \
 *           e:.5  e:.3  e:.2    e:.7  e:.2  e:.1      <- stage 2 (P(eps|s))
 *
 * writes it to a netCDF file, reads it back through the ScenarioGenerator
 * factory, and checks:
 *
 *  - the tree-navigation API (root / children / per-node stage, probability,
 *    data) — what an MSSB would use to build its Block tree;
 *  - the *view* API: several NodeView-s open at the same time iterate their
 *    own slice of the tree independently (no shared cursor) — the property
 *    that makes concurrent construction possible, and that lets a node-local
 *    consumer treat its slice as a plain ScenarioGenerator;
 *  - that joint leaf probabilities P(s)*P(eps|s) sum to 1 and match;
 *  - the legacy single-cursor walk still works (drop-in MSSG);
 *  - a serialize / deserialize round-trip.
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

#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using MSDSS = MultiStageDiscreteScenarioSet;

/*--------------------------------------------------------------------------*/

static int failures = 0;

static void check( bool cond , const std::string & what )
{
 std::cout << ( cond ? "  ok   " : "  FAIL " ) << what << std::endl;
 if( ! cond )
  ++failures;
 }

static void check_eq( double a , double b , const std::string & what )
{
 check( std::abs( a - b ) <= 1e-9 , what +
        " (got " + std::to_string( a ) + ", want " + std::to_string( b ) +
        ")" );
 }

/*--------------------------------------------------------------------------*/
// write the example tree to a netCDF file at the root group

static void write_tree( const std::string & fname )
{
 netCDF::NcFile f( fname , netCDF::NcFile::replace );

 f.putAtt( "type" , "MultiStageDiscreteScenarioSet" );

 auto sDim = f.addDim( "NumberStages" , 3 );
 auto nDim = f.addDim( "NumberNodes" , 9 );
 auto dDim = f.addDim( "ScenarioDataSize" , 1 );

 std::vector< unsigned int > ss = { 1 , 1 , 1 };
 f.addVar( "StageScenarioSize" , netCDF::NcUint() , sDim )
   .putVar( ss.data() );

 //                  root  s0    s1    e00   e01   e02   e10   e11   e12
 std::vector< unsigned int > stage  = { 0 , 1 , 1 , 2 , 2 , 2 , 2 , 2 , 2 };
 std::vector< unsigned int > parent = { 9 , 0 , 0 , 1 , 1 , 1 , 2 , 2 , 2 };
 std::vector< double > prob = { 1.0 , 0.6 , 0.4 ,
                               0.5 , 0.3 , 0.2 , 0.7 , 0.2 , 0.1 };
 std::vector< double > data = { 0 , 10 , 20 ,
                               100 , 101 , 102 , 200 , 201 , 202 };

 f.addVar( "NodeStage" , netCDF::NcUint() , nDim ).putVar( stage.data() );
 f.addVar( "NodeParent" , netCDF::NcUint() , nDim ).putVar( parent.data() );
 f.addVar( "NodeProbability" , netCDF::NcDouble() , nDim ).putVar(
                                                            prob.data() );
 f.addVar( "NodeData" , netCDF::NcDouble() , { nDim , dDim } ).putVar(
                                                            data.data() );
 }

/*--------------------------------------------------------------------------*/

int main( void )
{
 const std::string fname = "mssds_test.nc4";
 write_tree( fname );

 // read back through the ScenarioGenerator factory (exercises "type")
 auto * sg = ScenarioGenerator::deserialize( fname );
 check( sg != nullptr , "factory deserialize returns a generator" );
 auto * t = dynamic_cast< MSDSS * >( sg );
 check( t != nullptr , "it is a MultiStageDiscreteScenarioSet" );
 if( ! t ) { std::cerr << "fatal: cannot proceed" << std::endl; return( 1 ); }

 // -- tree navigation -----------------------------------------------------
 check( t->get_stage_number() == 3 , "3 stages" );
 const auto root = t->get_root();
 const auto & rc = t->get_children( root );
 check( rc.size() == 2 , "root has 2 children (s0, s1)" );
 check( t->get_node_stage( rc[ 0 ] ) == 1 , "children at stage 1" );
 check_eq( t->get_node_probability( rc[ 0 ] ) , 0.6 , "P(s0) = 0.6" );
 check_eq( t->get_node_probability( rc[ 1 ] ) , 0.4 , "P(s1) = 0.4" );

 // -- views in parallel (no shared cursor) --------------------------------
 // open one view per s-node at the same time and interleave their walks;
 // if a single cursor were shared they would clobber each other.
 auto * v0 = t->make_view( rc[ 0 ] );   // children of s0
 auto * v1 = t->make_view( rc[ 1 ] );   // children of s1

 check_eq( v0->get_current_scenario()[ 0 ] , 100 , "v0 child0 data" );
 check_eq( v1->get_current_scenario()[ 0 ] , 200 , "v1 child0 parallel" );
 check_eq( v0->get_current_scenario_probability() , 0.5 , "v0 P(eps|s0)" );
 check_eq( v1->get_current_scenario_probability() , 0.7 , "v1 P(eps|s1)" );

 check( v0->next_scenario() , "v0 advances" );
 check( v1->next_scenario() , "v1 advances" );
 check_eq( v0->get_current_scenario()[ 0 ] , 101 , "v0 child1 correct" );
 check_eq( v1->get_current_scenario()[ 0 ] , 201 , "v1 child1 correct" );

 check( v0->next_scenario() && ! v0->next_scenario() , "v0 has 3 children" );
 check( v1->next_scenario() && ! v1->next_scenario() , "v1 has 3 children" );
 delete v0;
 delete v1;

 // -- joint probabilities -------------------------------------------------
 double total = 0.0;
 double joint_s0_e0 = 0.0;
 for( auto s : t->get_children( root ) ) {
  const double ps = t->get_node_probability( s );
  for( auto e : t->get_children( s ) ) {
   const double j = ps * t->get_node_probability( e );
   total += j;
   if( ( s == t->get_children( root )[ 0 ] ) &&
       ( e == t->get_children( s )[ 0 ] ) )
    joint_s0_e0 = j;
   }
  }
 check_eq( total , 1.0 , "joint leaf probabilities sum to 1" );
 check_eq( joint_s0_e0 , 0.30 , "P(s0)*P(eps0|s0) = 0.30" );

 // -- legacy single cursor walk -------------------------------------------
 t->previous_stage( MultiStageScenarioGenerator::INFStage );
 check( t->get_current_stage() == 0 , "cursor back at stage 0" );
 check( t->next_stage() , "descend to stage 1" );
 check_eq( t->get_current_scenario()[ 0 ] , 10 , "cursor at s0" );
 check( t->next_scenario() , "move to sibling s1" );
 check_eq( t->get_current_scenario()[ 0 ] , 20 , "cursor at s1" );
 check( ! t->next_scenario() , "no more first-stage siblings" );
 check( t->next_stage() , "descend to stage 2 under s1" );
 check_eq( t->get_current_scenario()[ 0 ] , 200 , "cursor at first eps s1" );

 // -- round-trip ----------------------------------------------------------
 const std::string fname2 = "mssds_test_rt.nc4";
 {
  netCDF::NcFile f2( fname2 , netCDF::NcFile::replace );
  t->serialize( f2 );
 }
 auto * sg2 = ScenarioGenerator::deserialize( fname2 );
 auto * t2 = dynamic_cast< MSDSS * >( sg2 );
 check( t2 != nullptr , "round-trip: still a MultiStageDiscreteScenarioSet" );
 if( t2 ) {
  check( t2->get_stage_number() == 3 , "round-trip: 3 stages" );
  const auto r0 = t2->get_children( t2->get_root() )[ 0 ];
  check_eq( t2->get_node_probability( r0 ) , 0.6 ,
            "round-trip: P(s0) preserved" );
  }

 delete sg;
 delete sg2;
 std::remove( fname.c_str() );
 std::remove( fname2.c_str() );

 std::cout << "\n" << ( failures == 0 ? "All tests passed!!"
                                      : std::to_string( failures ) +
                                        " test(s) FAILED" ) << std::endl;
 return( failures == 0 ? 0 : 1 );
 }

/*--------------------------------------------------------------------------*/
/*----------------------- End File test_mssds.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
