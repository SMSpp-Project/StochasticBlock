/*--------------------------------------------------------------------------*/
/*-------------------- File DiscreteScenarioSet.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of DiscreteScenarioSet.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Benoît Tran \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Claude Opus 4.7 \n
 *         Antrophic \n
 *
 * \copyright &copy; by Antonio Frangioni, Benoît Tran
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DiscreteScenarioSet.h"
#include "Solver.h"

#include <cmath>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------- FACTORY REGISTRATION -------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( DiscreteScenarioSet );

/*--------------------------------------------------------------------------*/
/*------------------------ STATIC HELPER FUNCTIONS ------------------------*/
/*--------------------------------------------------------------------------*/
// Validate the poolSize parameter
static void validate_poolSize_value(
 DiscreteScenarioSet::ScenarioIndex size ,
 DiscreteScenarioSet::ScenarioIndex max_scenarios )
{
 if( size == 0 || size > max_scenarios )
  throw std::invalid_argument(
   "Invalid pool size: must be between 1 and " +
   std::to_string( max_scenarios ) );
 }

/*--------------------------------------------------------------------------*/

// Pick k distinct indices out of n, sampled with probability given by weights
static void generateWeightedRandomSubset(
 size_t n , size_t k ,
 const std::vector< double > & weights ,
 std::vector< ScenarioGenerator::ScenarioIndex > & ind ,
 std::mt19937 & rng )
{
 if( k > n )
  throw std::invalid_argument( "k must be <= n" );

 ind.clear();
 if( k == 0 ) return;

 if( k == n ) {
  ind.resize( n );
  std::iota( ind.begin() , ind.end() , 0 );
  return;
  }

 std::discrete_distribution< ScenarioGenerator::ScenarioIndex > dist(
  weights.begin() , weights.end() );

 std::unordered_set< ScenarioGenerator::ScenarioIndex > seen;
 seen.reserve( k );
 while( seen.size() < k )
  seen.insert( dist( rng ) );

 ind.assign( seen.begin() , seen.end() );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR ------------------------*/
/*--------------------------------------------------------------------------*/

DiscreteScenarioSet::DiscreteScenarioSet()
{
 set_seed( 1337 );
 }

/*--------------------------------------------------------------------------*/

DiscreteScenarioSet::~DiscreteScenarioSet()
{
 scenarioSet.resize( boost::extents[ 0 ][ 0 ] );
 scenarioIndexes.clear();
 scenarioIndexes.shrink_to_fit();
 universeIndexes.clear();
 universeIndexes.shrink_to_fit();
 setWeights.clear();
 setWeights.shrink_to_fit();
 poolWeights.clear();
 poolWeights.shrink_to_fit();
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ serialize --------------------------------*/
/*--------------------------------------------------------------------------*/

void DiscreteScenarioSet::serialize( netCDF::NcGroup & group ) const
{
 ScenarioGenerator::serialize( group );

 auto nbScenariosDim = group.addDim( "NumberScenarios" , nbScenarios );
 auto scenarioSizeDim = group.addDim( "ScenarioSize" , scenarioSize );

 auto scenariosVar = group.addVar(
  "Scenarios" , netCDF::NcDouble() ,
  { nbScenariosDim , scenarioSizeDim } );
 scenariosVar.putVar( scenarioSet.data() );

 auto probsVar = group.addVar(
  "PoolWeights" , netCDF::NcDouble() , nbScenariosDim );
 probsVar.putVar( setWeights.data() );
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ deserialize ------------------------------*/
/*--------------------------------------------------------------------------*/

void DiscreteScenarioSet::deserialize( const netCDF::NcGroup & group )
{
 currentScenarioIndex = 0;
 poolSize = 0;
 sumPoolWeights = 0.0;
 is_initialized = false;
 f_solution_ready = false;  // invalidate any cached solution

 scenarioSet.resize( boost::extents[ 0 ][ 0 ] );
 setWeights.clear();
 scenarioIndexes.clear();
 universeIndexes.clear();

 ::deserialize_dim( group , "NumberScenarios" , nbScenarios , false );
 ::deserialize_dim( group , "ScenarioSize" , scenarioSize , false );

 if( nbScenarios == 0 )
  throw std::invalid_argument( "NumberScenarios must be positive" );
 if( scenarioSize == 0 )
  throw std::invalid_argument( "ScenarioSize must be positive" );

 scenarioSet.resize( boost::extents[ nbScenarios ][ scenarioSize ] );
 std::vector< std::size_t > sizes = { nbScenarios , scenarioSize };
 ::deserialize( group , "Scenarios" , sizes , scenarioSet , true , false );

// Reserve space for setWeights before deserializing
 setWeights.reserve( nbScenarios );
// Deserialize probabilities (optional, default to uniform) 
 bool probsLoaded =
  ::deserialize( group , "PoolWeights" , nbScenarios , setWeights );

 if( ! probsLoaded )
  setWeights.assign( nbScenarios , 1.0 / nbScenarios );
 else if( setWeights.size() != nbScenarios )
  throw std::invalid_argument(
   "PoolWeights size does not match NumberScenarios" );

 double sum = std::accumulate( setWeights.begin() , setWeights.end() , 0.0 );
 if( std::abs( sum - 1.0 ) > 1e-6 )
  throw std::invalid_argument(
   "Scenario probabilities sum to " + std::to_string( sum ) +
   " rather than 1.0" );

 // Lazy canonical init: walk entire set in natural order
 universeIndexes.resize( nbScenarios );
 std::iota( universeIndexes.begin() , universeIndexes.end() , 0 );
 scenarioIndexes = universeIndexes;
 update_pool_weights();
 currentScenarioIndex = 0;
 is_initialized = true;
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------- load ----------------------------------*/
/*--------------------------------------------------------------------------*/

void DiscreteScenarioSet::load( std::istream & input )
{
 static const std::string sre( "DiscreteScenarioSet::load: stream read error" );

 // Read header: N scenarios of dimension D
 unsigned int N = 0 , D = 0;
 input >> N >> D;
 if( input.fail() )
  throw std::invalid_argument( sre + ": expected N D" );
 if( N == 0 ) throw std::invalid_argument( sre + ": N must be positive" );
 if( D == 0 ) throw std::invalid_argument( sre + ": D must be positive" );

 std::vector< double > weights( N );
 for( unsigned int i = 0 ; i < N ; ++i ) {
  input >> weights[ i ];
  if( input.fail() )
   throw std::invalid_argument( sre + ": failed reading weight " +
                                std::to_string( i ) );
  }

  double sum = std::accumulate( weights.begin() , weights.end() , 0.0 );
 if( std::abs( sum - 1.0 ) > 1e-6 )
  throw std::invalid_argument( sre + ": weights sum to " +
                               std::to_string( sum ) );

 std::vector< std::vector< double > > scenarios(
  N , std::vector< double >( D ) );
 for( unsigned int i = 0 ; i < N ; ++i )
  for( unsigned int d = 0 ; d < D ; ++d ) {
   input >> scenarios[ i ][ d ];
   if( input.fail() )
    throw std::invalid_argument( sre + ": failed reading scenario " +
                                 std::to_string( i ) );
   }

 load_from_memory( scenarios , weights );
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ load_from_memory -------------------------*/
/*--------------------------------------------------------------------------*/

void DiscreteScenarioSet::load_from_memory(
 const std::vector< std::vector< double > > & scenarios ,
 const std::vector< double > & weights )
{
 if( scenarios.empty() )
  throw std::invalid_argument(
   "DiscreteScenarioSet::load_from_memory: no scenarios provided" );

 const std::size_t N = scenarios.size();
 const std::size_t D = scenarios[ 0 ].size();
 for( std::size_t i = 1 ; i < N ; ++i )
  if( scenarios[ i ].size() != D )
   throw std::invalid_argument(
    "DiscreteScenarioSet::load_from_memory: scenario rows differ in size" );

 nbScenarios  = static_cast< ScenarioIndex >( N );
 scenarioSize = static_cast< ScenarioSize >( D );

 scenarioSet.resize( boost::extents[ N ][ D ] );
 for( std::size_t i = 0 ; i < N ; ++i )
  for( std::size_t d = 0 ; d < D ; ++d )
   scenarioSet[ i ][ d ] = scenarios[ i ][ d ];

 if( weights.empty() )
  setWeights.assign( N , 1.0 / static_cast< double >( N ) );
 else {
  if( weights.size() != N )
   throw std::invalid_argument(
    "DiscreteScenarioSet::load_from_memory: weights size mismatch" );
  setWeights = weights;
  }

 universeIndexes.clear();
 scenarioIndexes.clear();
 poolWeights.clear();
 sumPoolWeights = 0.0;
 currentScenarioIndex = 0;
 poolSize = 0;
 is_initialized = false;
 f_solution_ready = false;  // invalidate any cached solution
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ set_config --------------------------------*/
/*--------------------------------------------------------------------------*/

void DiscreteScenarioSet::set_config( Configuration * config )
{
 if( ! config ) return;

 // Pattern 1: SimpleConfiguration<int>, baseline, just pool size
 auto * simple = dynamic_cast< SimpleConfiguration< int > * >( config );
 if( simple ) {
  if( simple->f_value > 0 )
   poolSize = simple->f_value;
  return;
  }

 // Pattern 2: SimpleConfiguration<pair<int, Configuration*>> with solver
 auto * pair_cfg = dynamic_cast<
  SimpleConfiguration< std::pair< int , Configuration * > > * >( config );
 if( pair_cfg ) {
  if( pair_cfg->f_value.first > 0 )
   poolSize = pair_cfg->f_value.first;
  auto * sc = dynamic_cast< BlockSolverConfig * >(
   pair_cfg->f_value.second );
  if( sc )
   f_solver_config.reset( sc->clone() );
  return;
  }

 throw std::invalid_argument(
  "DiscreteScenarioSet::set_config: unsupported configuration type" );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- init_random_pool ----------------------------*/
/*--------------------------------------------------------------------------*/

// Build a pool by sampling pool_size scenarios at random, weighted by setWeights
void DiscreteScenarioSet::init_random_pool( ScenarioIndex pool_size )
{
 scenarioIndexes.clear();
 poolWeights.clear();
 sumPoolWeights = 0.0;
 currentScenarioIndex = 0;
 is_initialized = false;

 if( universeIndexes.empty() ) {
  universeIndexes.resize( nbScenarios );
  std::iota( universeIndexes.begin() , universeIndexes.end() , 0 );
  }

 const ScenarioIndex universe_size =
  static_cast< ScenarioIndex >( universeIndexes.size() );

 ScenarioIndex effective_size = pool_size;
 if( pool_size == INFScenario )
  effective_size = universe_size;
 else if( pool_size == 0 )
  throw std::invalid_argument(
   "DiscreteScenarioSet::init_random_pool: pool_size must be positive" );
 else if( pool_size > universe_size )
  throw std::invalid_argument(
   "DiscreteScenarioSet::init_random_pool: pool_size exceeds universe size" );

 poolSize = effective_size;

 if( effective_size > 0 ) {
  std::vector< double > universe_weights;
  universe_weights.reserve( universe_size );
  for( const auto idx : universeIndexes )
   universe_weights.push_back( setWeights[ idx ] );

  std::vector< ScenarioIndex > positions;
  generateWeightedRandomSubset(
   universe_size , effective_size , universe_weights , positions , rng );

  scenarioIndexes.reserve( effective_size );
  for( const auto p : positions )
   scenarioIndexes.push_back( universeIndexes[ p ] );

  update_pool_weights();
  is_initialized = true;
  }
 }

/*--------------------------------------------------------------------------*/
/*------------------------ init_representative_pool -----------------------*/
/*--------------------------------------------------------------------------*/

void DiscreteScenarioSet::init_representative_pool(
 ScenarioIndex target_pool_size )
{
 scenarioIndexes.clear();
 poolWeights.clear();
 sumPoolWeights = 0.0;
 currentScenarioIndex = 0;
 is_initialized = false;

 // INFScenario means no restriction, reset to full set
 if( target_pool_size == INFScenario ||
     target_pool_size == nbScenarios ) {
  universeIndexes.resize( nbScenarios );
  std::iota( universeIndexes.begin() , universeIndexes.end() , 0 );
  scenarioIndexes = universeIndexes;
  poolSize = nbScenarios;
  update_pool_weights();
  is_initialized = ! scenarioIndexes.empty();
  return;
  }

 if( target_pool_size == 0 )
  throw std::invalid_argument(
   "DiscreteScenarioSet::init_representative_pool: "
   "target_pool_size must be positive" );

 if( target_pool_size > nbScenarios )
  throw std::invalid_argument(
   "DiscreteScenarioSet::init_representative_pool: "
   "target_pool_size exceeds nbScenarios" );

 poolSize = target_pool_size;

 if( f_solution_ready ) {
  currentScenarioIndex = 0;
  universeIndexes = scenarioIndexes;
  if( ! scenarioIndexes.empty() )
   is_initialized = true;
  return;
  }

 if( f_solver_config ) {
  // Create ScenarioReductionBlock
  auto srb = std::make_unique< ScenarioReductionBlock >();
  srb->set_scenario_generator( this );

  if( f_stochastic_block )
   srb->set_stochastic_block( f_stochastic_block );

  // Attach Solver via BlockSolverConfig
  auto config_copy = std::unique_ptr< BlockSolverConfig >(
   f_solver_config->clone() );
  config_copy->apply( srb.get() );

  if( srb->get_registered_solvers().empty() )
   throw std::runtime_error(
    "DiscreteScenarioSet::init_representative_pool: "
    "no Solver registered after BlockSolverConfig::apply" );

  auto * solver = srb->get_registered_solvers().front();

  // Solve
  int status = solver->compute();
  if( status != Solver::kOK )
   throw std::runtime_error(
    "DiscreteScenarioSet::init_representative_pool: "
    "Solver failed with status " + std::to_string( status ) );

  solver->get_var_solution();

  // Read solution from ScenarioReductionBlock
  const auto & sol = srb->get_solution();
  if( ! sol.is_set() )
   throw std::runtime_error(
    "DiscreteScenarioSet::init_representative_pool: "
    "Solver did not write a solution into ScenarioReductionBlock" );

  scenarioIndexes = sol.selected_indices;
  update_pool_weights_with_assignments( sol.assignments );
  f_solution_ready = true;  // cache solution for subsequent calls
  }
 else {
  // No solver, fall back to baseline selection
  apply_baseline_selection( target_pool_size );
  update_pool_weights();
  }

 poolSize = scenarioIndexes.size();
 universeIndexes = scenarioIndexes;
 currentScenarioIndex = 0;
 if( ! scenarioIndexes.empty() )
  is_initialized = true;
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ next_scenario ----------------------------*/
/*--------------------------------------------------------------------------*/

// Advance to the next scenario in the pool; false when already at the last one.
bool DiscreteScenarioSet::next_scenario()
{
 if( ! is_initialized || scenarioIndexes.empty() )
  return false;

 if( currentScenarioIndex < scenarioIndexes.size() - 1 ) {
  ++currentScenarioIndex;
  return true;
  }

 return false;
 }

/*--------------------------------------------------------------------------*/
/*------------------------ get_selected_scenarios -------------------------*/
/*--------------------------------------------------------------------------*/

std::span< const DiscreteScenarioSet::ScenarioIndex >
DiscreteScenarioSet::get_selected_scenarios() const
{
 if( ! is_initialized )
  throw std::runtime_error(
   "get_selected_scenarios: pool not initialized" );
 return scenarioIndexes;
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- HELPER METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void DiscreteScenarioSet::empty_pool()
{
 scenarioIndexes.clear();
 scenarioIndexes.shrink_to_fit();
 poolWeights.clear();
 poolWeights.shrink_to_fit();
 currentScenarioIndex = 0;
 sumPoolWeights = 0.0;
 is_initialized = false;
 }

/*--------------------------------------------------------------------------*/

// Symmetric matrix of pairwise distances between all scenarios (the transport costs)
DiscreteScenarioSet::CMatrix
DiscreteScenarioSet::compute_transport_cost_matrix(
 ScenarioIndex n_scenarios , ScenarioSize scenario_size , float ell ) const
{
 CMatrix costs( boost::extents[ n_scenarios ][ n_scenarios ] );

 Eigen::Map< const Eigen::Matrix< double , Eigen::Dynamic ,
  Eigen::Dynamic , Eigen::RowMajor > >
  all_scenarios( scenarioSet.data() , n_scenarios , scenario_size );

 for( ScenarioIndex i = 0 ; i < n_scenarios ; ++i ) {
  costs[ i ][ i ] = 0.0;
  for( ScenarioIndex j = i + 1 ; j < n_scenarios ; ++j ) {
   double d = compute_scenario_distance(
    all_scenarios.row( i ) , all_scenarios.row( j ) , ell );
   costs[ i ][ j ] = d;
   costs[ j ][ i ] = d;
   }
  }

 return costs;
 }

/*--------------------------------------------------------------------------*/

double DiscreteScenarioSet::compute_scenario_distance(
 const Eigen::VectorXd & s1 , const Eigen::VectorXd & s2 , float ell ) const
{
 return std::pow( ( s1 - s2 ).norm() , ell );
 }

/*--------------------------------------------------------------------------*/

// Renormalize the weights of the pooled scenarios so they sum to 1
void DiscreteScenarioSet::update_pool_weights()
{
 sumPoolWeights = 0.0;
 for( const auto & idx : scenarioIndexes )
  sumPoolWeights += setWeights[ idx ];

 poolWeights.clear();
 poolWeights.reserve( scenarioIndexes.size() );

 if( sumPoolWeights > 0.0 )
  for( const auto & idx : scenarioIndexes )
   poolWeights.push_back( setWeights[ idx ] / sumPoolWeights );
 else {
  double u = 1.0 / scenarioIndexes.size();
  poolWeights.assign( scenarioIndexes.size() , u );
  }
 }

/*--------------------------------------------------------------------------*/

// Each representative gets the total weight of the scenarios assigned to it,
// then the pool weights are renormalized to sum to 1
void DiscreteScenarioSet::update_pool_weights_with_assignments(
 const std::vector< ScenarioIndex > & assignments )
{
 if( assignments.size() != nbScenarios )
  throw std::invalid_argument(
   "update_pool_weights_with_assignments: assignments size mismatch" );

 std::unordered_map< ScenarioIndex , size_t > repr_to_pos;
 for( size_t pos = 0 ; pos < scenarioIndexes.size() ; ++pos )
  repr_to_pos[ scenarioIndexes[ pos ] ] = pos;

 std::vector< double > agg( scenarioIndexes.size() , 0.0 );
 for( ScenarioIndex i = 0 ; i < nbScenarios ; ++i ) {
  auto it = repr_to_pos.find( assignments[ i ] );
  if( it != repr_to_pos.end() )
   agg[ it->second ] += setWeights[ i ];
  else
   throw std::runtime_error(
    "update_pool_weights_with_assignments: scenario " +
    std::to_string( i ) + " assigned to non-selected representative" );
  }

 sumPoolWeights = std::accumulate( agg.begin() , agg.end() , 0.0 );

 poolWeights.clear();
 poolWeights.reserve( scenarioIndexes.size() );

 if( sumPoolWeights > 0.0 )
  for( const auto & w : agg )
   poolWeights.push_back( w / sumPoolWeights );
 else {
  double u = 1.0 / scenarioIndexes.size();
  poolWeights.assign( scenarioIndexes.size() , u );
  }
 }

/*--------------------------------------------------------------------------*/

// Fallback selection without a Solver: just keep the target_size heaviest scenarios
void DiscreteScenarioSet::apply_baseline_selection( ScenarioIndex target_size )
{
 std::vector< std::pair< ScenarioIndex , double > > wp;
 wp.reserve( nbScenarios );
 for( ScenarioIndex i = 0 ; i < nbScenarios ; ++i )
  wp.emplace_back( i , ( i < setWeights.size() ?
   setWeights[ i ] : 1.0 / nbScenarios ) );

 std::sort( wp.begin() , wp.end() ,
  []( const auto & a , const auto & b ) { return a.second > b.second; } );

 scenarioIndexes.clear();
 scenarioIndexes.reserve( target_size );
 for( ScenarioIndex i = 0 ; i < target_size && i < nbScenarios ; ++i )
  scenarioIndexes.push_back( wp[ i ].first );

 std::sort( scenarioIndexes.begin() , scenarioIndexes.end() );
 }

/*--------------------------------------------------------------------------*/
/*-------------- End file DiscreteScenarioSet.cpp --------------------------*/
/*--------------------------------------------------------------------------*/