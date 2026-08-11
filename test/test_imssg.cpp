/*--------------------------------------------------------------------------*/
/*-------------------------- File test_imssg.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Stand-alone test / demo for IndependentMultiStageScenarioGenerator, the
 * stage-independent MultiStageScenarioGenerator, and specifically for its
 * View: the same view-based interface that a genuine scenario *tree*
 * implements has to make sense for the degenerate case where the stages do
 * not depend on the history at all, which is what this test checks.
 *
 * It builds a 3-stage generator whose stages hold 2, 3 and 2 realizations
 * respectively (one DiscreteScenarioSet each), writes it to a netCDF file,
 * reads it back through the ScenarioGenerator factory, and checks:
 *
 *  - that the tree is walked with the same View API as any other
 *    :MultiStageScenarioGenerator: root_view() + descend() + climb(), each
 *    of the last two undoing the other;
 *  - that the pool a View reads is the one of the stage it has reached, and
 *    that, the stages being independent, descend() does not depend on which
 *    realization is currently selected;
 *  - that views pinned at *distinct* stages are independent of one another,
 *    while two views pinned at the *same* stage share the cursor held by
 *    that stage's inner :ScenarioGenerator (see the class comments);
 *  - that the generator, seen as a plain ScenarioGenerator, is its own root
 *    view.
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

#include "IndependentMultiStageScenarioGenerator.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using IMSSG = IndependentMultiStageScenarioGenerator;

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
// write one DiscreteScenarioSet, with one-component scenarios, in a group

static void write_stage( netCDF::NcGroup & g ,
                         const std::vector< double > & values )
{
 g.putAtt( "type" , "DiscreteScenarioSet" );

 const auto n = values.size();
 auto nDim = g.addDim( "NumberScenarios" , n );
 auto sDim = g.addDim( "ScenarioSize" , 1 );

 g.addVar( "Scenarios" , netCDF::NcDouble() , { nDim , sDim } ).putVar(
                                                          values.data() );

 std::vector< double > w( n , 1.0 / double( n ) );
 g.addVar( "PoolWeights" , netCDF::NcDouble() , nDim ).putVar( w.data() );
 }

/*--------------------------------------------------------------------------*/
// write a 3-stage independent generator: 2, 3 and 2 realizations

static void write_generator( const std::string & fname )
{
 netCDF::NcFile f( fname , netCDF::NcFile::replace );

 f.putAtt( "type" , "IndependentMultiStageScenarioGenerator" );
 f.addDim( "NumberStages" , 3 );

 auto g0 = f.addGroup( "Stage_0" );
 write_stage( g0 , { 10 , 20 } );
 auto g1 = f.addGroup( "Stage_1" );
 write_stage( g1 , { 100 , 200 , 300 } );
 auto g2 = f.addGroup( "Stage_2" );
 write_stage( g2 , { 1000 , 2000 } );
 }

/*--------------------------------------------------------------------------*/

int main( void )
{
 const std::string fname = "imssg_test.nc4";
 write_generator( fname );

 // read back through the ScenarioGenerator factory (exercises "type")
 auto * sg = ScenarioGenerator::deserialize( fname );
 check( sg != nullptr , "factory deserialize returns a generator" );
 auto * g = dynamic_cast< IMSSG * >( sg );
 check( g != nullptr , "it is an IndependentMultiStageScenarioGenerator" );
 if( ! g ) { std::cerr << "fatal: cannot proceed" << std::endl; return( 1 ); }

 check( g->get_stage_number() == 3 , "3 stages" );
 check( g->is_stage_independent() , "the stages are independent" );

 // -- the root view -------------------------------------------------------
 auto root = g->root_view();
 check( root->stage() == 0 , "the root view is pinned at stage 0" );
 check( root->get_support_size() == 2 , "stage 0 has 2 realizations" );
 check_eq( root->get_current_scenario()[ 0 ] , 10 , "the first of them" );
 check_eq( root->get_current_scenario_probability() , 0.5 , "P = 1/2" );
 check( root->next_scenario() , "the root view advances" );
 check_eq( root->get_current_scenario()[ 0 ] , 20 , "the second one" );
 check( ! root->next_scenario() , "no third realization at stage 0" );

 // -- descend / climb -----------------------------------------------------
 auto w = root->clone();
 check( w->descend() , "descend to stage 1" );
 check( w->stage() == 1 , "now pinned at stage 1" );
 check( w->get_support_size() == 3 , "stage 1 has 3 realizations" );
 check_eq( w->get_current_scenario()[ 0 ] , 100 ,
           "descend() selects the first realization of the new stage" );
 check( w->descend() , "descend to stage 2" );
 check( w->get_support_size() == 2 , "stage 2 has 2 realizations" );
 check_eq( w->get_current_scenario()[ 0 ] , 1000 , "the first of them" );
 check( ! w->descend() , "stage 2 is the last one" );
 check( w->climb() , "climb back to stage 1" );
 check( w->stage() == 1 , "back at stage 1" );
 check( w->climb() && w->stage() == 0 , "climb back to stage 0" );
 check( ! w->climb() , "cannot climb past the first stage" );

 // -- the history is immaterial -------------------------------------------
 // the stages being independent, descending after having selected another
 // realization must land on the very same pool
 auto u = g->root_view();
 check( u->next_scenario() , "select the other realization of stage 0" );
 check( u->descend() , "descend from it" );
 check( u->get_support_size() == 3 , "the stage-1 pool is the same" );
 check_eq( u->get_current_scenario()[ 0 ] , 100 , "and so is its first" );

 // -- views at distinct stages are independent ----------------------------
 auto v1 = g->root_view();
 check( v1->descend() , "v1 goes to stage 1" );
 auto v2 = g->root_view();
 check( v2->descend() && v2->descend() , "v2 goes to stage 2" );

 check( v1->next_scenario() , "v1 advances at stage 1" );
 check_eq( v1->get_current_scenario()[ 0 ] , 200 , "v1 second of stage 1" );
 check_eq( v2->get_current_scenario()[ 0 ] , 1000 ,
           "v2 is unaffected: stages are separate objects" );
 check( v2->next_scenario() , "v2 advances at stage 2" );
 check_eq( v2->get_current_scenario()[ 0 ] , 2000 , "v2 second of stage 2" );
 check_eq( v1->get_current_scenario()[ 0 ] , 200 , "v1 still where it was" );

 // two views pinned at the SAME stage share that stage's inner cursor:
 // this is a property of this generator, whose inners hold the cursor,
 // not of the View interface
 auto v3 = v1->clone();
 check_eq( v3->get_current_scenario()[ 0 ] , 200 ,
           "a clone starts where the view it was taken from is" );
 check( v3->next_scenario() , "the clone advances at stage 1" );
 check_eq( v1->get_current_scenario()[ 0 ] , 300 ,
           "and the other view at the same stage follows: shared cursor" );

 // -- the generator is its own root view ----------------------------------
 v1->climb();
 auto r = g->root_view();
 r->reset_pool();
 check( g->get_support_size() == 2 , "the generator reads stage 0" );
 check_eq( g->get_current_scenario()[ 0 ] , 10 , "the generator is at 10" );
 check( g->next_scenario() , "the generator advances" );
 check_eq( g->get_current_scenario()[ 0 ] , 20 , "the generator is at 20" );
 g->reset_pool();
 check_eq( g->get_current_scenario()[ 0 ] , 10 , "the generator rewinds" );

 // -- a pool re-definition invalidates the views that relied on it --------
 // the stages being independent, re-defining the pool of one of them only
 // concerns the views pinned there, and not the one that re-defined it
 auto a = g->root_view();
 auto b = g->root_view();
 auto c = g->root_view();
 check( c->descend() , "c moves to stage 1" );
 check( a->is_valid() && b->is_valid() && c->is_valid() ,
        "all the views are valid to start with" );

 a->init_random_pool();                   // re-defines the pool of stage 0
 check( a->is_valid() , "the view that re-defined the pool is still valid" );
 check( ! b->is_valid() , "the other view at that stage is invalidated" );
 check( c->is_valid() , "a view at another stage is not" );

 bool refused = false;
 try { b->get_support_size(); } catch( const std::exception & ) {
  refused = true;
  }
 check( refused , "an invalidated view refuses to be read" );

 refused = false;
 try { b->descend(); } catch( const std::exception & ) { refused = true; }
 check( refused , "and refuses to be moved" );

 delete sg;
 std::remove( fname.c_str() );

 std::cout << "\n" << ( failures == 0 ? "All tests passed!!"
                                      : std::to_string( failures ) +
                                        " test(s) FAILED" ) << std::endl;
 return( failures == 0 ? 0 : 1 );
 }

/*--------------------------------------------------------------------------*/
/*----------------------- End File test_imssg.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
