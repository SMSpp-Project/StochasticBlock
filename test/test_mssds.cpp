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
 *  - the view API, the only way the tree is read: the pool of a View are the
 *    realizations available at the position it pins, and descend() / climb()
 *    move it through the tree, each undoing the other;
 *  - that a clone() moves independently of the View it was taken from, and
 *    that several views open at the same time iterate their own slice of the
 *    tree without interfering — the property that makes concurrent
 *    construction possible, and that lets a node-local consumer treat its
 *    slice as a plain ScenarioGenerator;
 *  - that joint leaf probabilities P(s)*P(eps|s) sum to 1 and match;
 *  - that the generator, seen as a plain ScenarioGenerator, is its own root
 *    view;
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

 // -- the root view -------------------------------------------------------
 check( t->get_stage_number() == 3 , "3 stages" );
 auto root = t->root_view();
 check( root->stage() == 0 , "the root view is pinned at stage 0" );
 check( root->get_support_size() == 2 , "its pool are s0 and s1" );
 check_eq( root->get_current_scenario()[ 0 ] , 10 , "first realization s0" );
 check_eq( root->get_current_scenario_probability() , 0.6 , "P(s0) = 0.6" );
 check( root->next_scenario() , "the root view advances to s1" );
 check_eq( root->get_current_scenario()[ 0 ] , 20 , "second one is s1" );
 check_eq( root->get_current_scenario_probability() , 0.4 , "P(s1) = 0.4" );
 check( ! root->next_scenario() , "no third first-stage realization" );
 root->reset_pool();
 check_eq( root->get_current_scenario()[ 0 ] , 10 , "reset_pool() rewinds" );

 // -- descend / climb -----------------------------------------------------
 // a clone moves on its own; descend() extends the pinned history with the
 // realization currently selected, climb() drops it again
 auto w = root->clone();
 check( w->next_scenario() , "the clone moves to s1" );
 check_eq( root->get_current_scenario()[ 0 ] , 10 ,
           "the clone moved, the view it was taken from did not" );
 check( w->descend() , "the clone descends into s1" );
 check( w->stage() == 1 , "it is now pinned at stage 1" );
 check( w->get_support_size() == 3 , "s1 has 3 inner realizations" );
 check_eq( w->get_current_scenario()[ 0 ] , 200 , "the first of them" );
 check_eq( w->get_current_scenario_probability() , 0.7 , "P(eps0|s1)" );
 check( ! w->clone()->descend() , "the stage-2 realizations are leaves" );
 check( w->climb() , "the clone climbs back" );
 check( w->stage() == 0 , "back at stage 0" );
 check_eq( w->get_current_scenario()[ 0 ] , 20 ,
           "climb() re-selects the realization descend() went into" );
 check( ! w->climb() , "cannot climb past the root" );

 // -- views used at the same time (no shared cursor) ----------------------
 // open one view per s-node at the same time and interleave their walks;
 // if a single cursor were shared they would clobber each other
 auto v0 = root->clone();
 check( v0->descend() , "v0 descends into s0" );
 auto v1 = root->clone();
 check( v1->next_scenario() , "v1 selects s1" );
 check( v1->descend() , "v1 descends into s1" );

 check_eq( v0->get_current_scenario()[ 0 ] , 100 , "v0 realization 0" );
 check_eq( v1->get_current_scenario()[ 0 ] , 200 , "v1 realization 0" );
 check_eq( v0->get_current_scenario_probability() , 0.5 , "v0 P(eps|s0)" );
 check_eq( v1->get_current_scenario_probability() , 0.7 , "v1 P(eps|s1)" );

 check( v0->next_scenario() , "v0 advances" );
 check( v1->next_scenario() , "v1 advances" );
 check_eq( v0->get_current_scenario()[ 0 ] , 101 , "v0 realization 1" );
 check_eq( v1->get_current_scenario()[ 0 ] , 201 , "v1 realization 1" );

 check( v0->next_scenario() && ! v0->next_scenario() , "v0 has 3 of them" );
 check( v1->next_scenario() && ! v1->next_scenario() , "v1 has 3 of them" );

 // -- joint probabilities -------------------------------------------------
 double total = 0.0;
 double joint_s0_e0 = 0.0;
 auto s = t->root_view();
 bool first_s = true;
 do {
  const double ps = s->get_current_scenario_probability();
  auto e = s->clone();
  if( ! e->descend() )
   continue;
  bool first_e = true;
  do {
   const double j = ps * e->get_current_scenario_probability();
   total += j;
   if( first_s && first_e )
    joint_s0_e0 = j;
   first_e = false;
   } while( e->next_scenario() );
  first_s = false;
  } while( s->next_scenario() );

 check_eq( total , 1.0 , "joint leaf probabilities sum to 1" );
 check_eq( joint_s0_e0 , 0.30 , "P(s0)*P(eps0|s0) = 0.30" );

 // -- the generator is its own root view ----------------------------------
 check( t->get_support_size() == 2 , "the generator reads the root pool" );
 check_eq( t->get_current_scenario()[ 0 ] , 10 , "the generator is at s0" );
 check( t->next_scenario() , "the generator moves to s1" );
 check_eq( t->get_current_scenario()[ 0 ] , 20 , "the generator is at s1" );
 check( ! t->next_scenario() , "no third first-stage realization" );
 t->reset_pool();
 check_eq( t->get_current_scenario()[ 0 ] , 10 , "the generator rewinds" );

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
  auto r2 = t2->root_view();
  check( r2->get_support_size() == 2 , "round-trip: root still has 2" );
  check_eq( r2->get_current_scenario_probability() , 0.6 ,
            "round-trip: P(s0) preserved" );
  check( r2->descend() , "round-trip: still walkable" );
  check_eq( r2->get_current_scenario()[ 0 ] , 100 ,
            "round-trip: inner realization preserved" );
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
