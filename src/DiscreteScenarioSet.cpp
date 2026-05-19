/*--------------------------------------------------------------------------*/
/*-------------------- File DiscreteScenarioSet.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Class DiscreteScenarioSet that is an implementation of ScenarioGenerator
 * suited to the case where the input distribution is contained in a netCDF
 * file as a collection of vectors. Associated with the header file
 * DiscreteScenarioSet.h
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Benoit Tran \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Claude Opus 4.7 \n
 *         Antrophic \n
 *
 * \copyright &copy; by Antonio Frangioni and Benoit Tran
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DiscreteScenarioSet.h"

#include "CapacitatedFacilityLocationBlock.h"
#include "Solver.h"

#include <chrono>         // For timestamp generation
#include <cmath>          // For std::abs
#include <cstdio>         // For std::remove
#include <numeric>        // For std::accumulate (not in SMSTypedefs)
#include <unordered_set>  // For std::unordered_set (rejection sampling)
#include <unordered_map>  // For std::unordered_map (weight aggregation)

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------- FACTORY MANAGEMENT ----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( DiscreteScenarioSet );

/*--------------------------------------------------------------------------*/
/*------------------------- STATIC HELPER FUNCTIONS ------------------------*/
/*--------------------------------------------------------------------------*/
/* Static Helper Functions
 * These are internal helper functions not part of the public interface.
 * They provide utility functionality used by the class implementation. */

// Validate the poolSize parameter
static void
validate_poolSize_value( DiscreteScenarioSet::ScenarioIndex size,
                         DiscreteScenarioSet::ScenarioIndex max_scenarios )
{
 if( size == 0 || size > max_scenarios ) {
  throw std::invalid_argument(
      "Invalid pool size parameter: must be between 1 and " +
      std::to_string( max_scenarios ) );
 }
}

/*--------------------------------------------------------------------------*/
// Draw k elements among n with weighted sampling
/* The function generateWeightedRandomSubset draws k elements among n using
 * weighted random sampling without replacement. The chosen indexes are
 * moved into the input variable ind. */

static void generateWeightedRandomSubset( size_t n , size_t k ,
				       const std::vector<double> & weights ,
		       std::vector<ScenarioGenerator::ScenarioIndex> & ind ,
					  std::mt19937 & rng )
{
 if( k > n )
  throw( std::invalid_argument( "k must be less or equal than n." ) );

 // Clear the output container
 ind.clear();

 // Special case: if k is 0, just return empty vector
 if( k == 0 ) {
  return;
 }

 // Special case: if k == n, return all indices
 if( k == n ) {
  ind.resize( n );
  std::iota( ind.begin(), ind.end(), 0 );
  return;
 }

 // Use weighted sampling without replacement via rejection sampling
 std::discrete_distribution<ScenarioGenerator::ScenarioIndex> dist(
     weights.begin(), weights.end() );

 // Use unordered_set to track unique selections
 std::unordered_set<ScenarioGenerator::ScenarioIndex> selected_set;
 selected_set.reserve( k );

 // Keep sampling until we have k unique indices
 while ( selected_set.size() < k ) {
  ScenarioGenerator::ScenarioIndex selected = dist( rng );
  selected_set.insert( selected );  // Automatically handles duplicates
 }

 // Convert set to vector
 ind.clear();
 ind.reserve( k );
 ind.assign( selected_set.begin(), selected_set.end() );
 }

/*--------------------------------------------------------------------------*/
// Generate default BlockConfig for CFL (static helper function)
static BlockConfig * generate_default_cfl_config(
				    DiscreteScenarioSet::ScenarioIndex size )
{
 auto * config = new BlockConfig( false );  // not differential

 // Set pool size in extra_Configuration
 config->f_extra_Configuration = new SimpleConfiguration<int>( size );

 return config;
 }

/*--------------------------------------------------------------------------*/
// Create CFL problem data structures
static std::tuple< CapacitatedFacilityLocationBlock::DVector ,
                   CapacitatedFacilityLocationBlock::CVector ,
                   CapacitatedFacilityLocationBlock::DVector >
create_cfl_problem_data( DiscreteScenarioSet::ScenarioIndex n_scenarios ,
                         const std::vector<double> & setWeights )
{
 CapacitatedFacilityLocationBlock::DVector capacities( n_scenarios );
 CapacitatedFacilityLocationBlock::CVector fixed_costs( n_scenarios );
 CapacitatedFacilityLocationBlock::DVector demands( n_scenarios );

 // Set up the basic CFL parameters
 for( DiscreteScenarioSet::ScenarioIndex i = 0; i < n_scenarios; ++i ) {
  capacities[i] = 1.0;   // Capacities = weights of the reduced distribution
  fixed_costs[i] = 0.0;  // No fixed cost in scenario reduction
  demands[i] = setWeights[i];  // Demand equals weight
  }

 return std::make_tuple( std::move( capacities ), std::move( fixed_costs ),
                         std::move( demands ) );
 }

/*--------------------------------------------------------------------------*/
// Extract selected scenarios and assignments from CFL block results
static void extract_scenarios_and_assignments_from_cfl_block(
    const CapacitatedFacilityLocationBlock * cflBlock ,
    DiscreteScenarioSet::ScenarioIndex n_scenarios ,
    std::vector<DiscreteScenarioSet::ScenarioIndex> & scenarioIndexes ,
    std::vector<DiscreteScenarioSet::ScenarioIndex> & scenario_assignments )
{
#ifndef NDEBUG
 std::cout << "DEBUG [extract_scenarios_and_assignments_from_cfl_block]: "
              "Starting extraction from CFL block"
           << std::endl;
#endif

 // Clear existing selection and reserve space
 scenarioIndexes.clear();
 scenarioIndexes.reserve( n_scenarios );

 // Initialize assignments vector
 scenario_assignments.clear();
 scenario_assignments.resize( n_scenarios );

 // identify open facilities (selected scenarios)
 std::vector<bool> is_open( n_scenarios, false );
 for( DiscreteScenarioSet::ScenarioIndex i = 0; i < n_scenarios; ++i ) {
  // Get the y variable for facility i
  const auto * y_var = cflBlock->get_y( i );
  if( !y_var ) {
   throw std::runtime_error( "Failed to get y variable for facility " +
                             std::to_string( i ) );
  }

  // Get value from the variable directly
  double y_value = y_var->get_value();
#ifndef NDEBUG
  std::cout << "DEBUG [extract_scenarios_and_assignments_from_cfl_block]: y["
            << i << "] = " << y_value << std::endl;
#endif

  // Check if facility i is open (y[i] > 0.5)
  if( y_value > 0.5 ) {
   scenarioIndexes.push_back( i );
   is_open[i] = true;
  }
 }

 // determine assignments (which facility serves each customer)
 // In scenario reduction with unsplittable demand, each scenario is assigned to
 // exactly one representative
 for( DiscreteScenarioSet::ScenarioIndex customer = 0; customer < n_scenarios;
       ++customer ) {
  DiscreteScenarioSet::ScenarioIndex assigned_facility = 0;
  bool found_assignment = false;

  // Find the single facility this customer is assigned to
  // (x[facility][customer] = 1.0)
  for( DiscreteScenarioSet::ScenarioIndex facility = 0; facility < n_scenarios;
        ++facility ) {
   if( !is_open[facility] )
    continue;

   // Get the x variable for assignment of customer to facility
   const auto * x_var = cflBlock->get_x( facility, customer );
   if( !x_var )
    continue;

   double x_value = x_var->get_value();

   // With unsplittable demand, x should be either 0 or 1
   if( x_value >
        0.5 ) {  // Consider values > 0.5 as 1 (for numerical tolerance)
    assigned_facility = facility;
    found_assignment = true;
    break;  // Found the single assignment, no need to check others
   }
  }

  if( !found_assignment ) {
   // This shouldn't happen if the solver correctly sets x variables
   throw std::runtime_error(
       "No assignment found for scenario " + std::to_string( customer ) +
       ". The solver may not have set x variables correctly." );
  }

  scenario_assignments[customer] = assigned_facility;

#ifndef NDEBUG
  if( found_assignment &&
       customer < 10 ) {  // Only print first 10 for debugging
   std::cout
       << "DEBUG [extract_scenarios_and_assignments_from_cfl_block]: Scenario "
       << customer << " assigned to representative " << assigned_facility
       << std::endl;
  }
#endif
 }

#ifndef NDEBUG
 std::cout
     << "DEBUG [extract_scenarios_and_assignments_from_cfl_block]: Selected "
     << scenarioIndexes.size() << " representative scenarios: ";
 for( auto idx : scenarioIndexes )
  std::cout << idx << " ";
 std::cout << std::endl;
#endif

 // If no scenarios were selected, throw an error
 if( scenarioIndexes.empty() ) {
  throw std::runtime_error(
      "No scenarios selected by the reduction algorithm" );
 }
}

/*--------------------------------------------------------------------------*/
/*----------------------- GETTERS AND SETTERS ------------------------------*/
/*--------------------------------------------------------------------------*/

// Get the indices of currently selected scenarios
std::span< const DiscreteScenarioSet::ScenarioIndex >
DiscreteScenarioSet::get_selected_scenarios( void ) const
{
 if( ! is_initialized )
  throw( std::runtime_error( "Pool has not been initialized: call "
			     "init_random_pool() or "
			     "init_representative_pool() first" ) );

 return scenarioIndexes;  // Implicit conversion to span
 }

/*--------------------------------------------------------------------------*/
/*-------------------- HELPER METHODS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

void DiscreteScenarioSet::empty_pool( void )
{
 // Clear the indexes and free memory
 scenarioIndexes.clear();
 scenarioIndexes.shrink_to_fit();

 // Clear pool weights and free memory
 poolWeights.clear();
 poolWeights.shrink_to_fit();

 // Reset pool state (poolSize is configuration, not actual size)
 currentScenarioIndex = 0;
 // poolSize is NOT reset - it's the configured/desired size
 sumPoolWeights = 0.0;
 is_initialized = false;
 }

/*--------------------------------------------------------------------------*/
/*------ SCENARIO REDUCTION CONFIGURATION METHODS --------------------------*/
/*--------------------------------------------------------------------------*/

// Set scenario reduction configuration using mixed ownership model
void DiscreteScenarioSet::set_config( BlockConfig * block_config,
                                      BlockSolverConfig * solver_config )
{
 // Clean up existing configurations
 if( f_block_config )
  delete f_block_config;

 // f_solver_config will be cleaned up automatically by unique_ptr

 // Set the new configurations following SMS++ ownership patterns:
 // - Clone the BlockConfig (SMS++ pattern: blocks clone their configs)
 // - Take ownership of BlockSolverConfig via unique_ptr
 f_block_config = block_config ? block_config->clone() : nullptr;
 f_solver_config.reset( solver_config );

 // Extract and set poolSize from the BlockConfig if available
 if( block_config && block_config->f_extra_Configuration ) {
  auto * poolSize_config = dynamic_cast<SimpleConfiguration<int> *>(
      block_config->f_extra_Configuration );
  if( poolSize_config && poolSize_config->f_value > 0 )
   poolSize = poolSize_config->f_value;
  }
 }

/*--------------------------------------------------------------------------*/
// Set scenario reduction configuration with poolSize parameter

void DiscreteScenarioSet::set_config( BlockConfig * block_config,
                                      BlockSolverConfig * solver_config,
                                      ScenarioIndex pool_size )
{
 // Validate and set pool size
 validate_poolSize_value( pool_size, nbScenarios );
 poolSize = pool_size;

 // Then set the configurations
 set_config( block_config, solver_config );
 }

/*--------------------------------------------------------------------------*/
// Set configuration from a Configuration object

void DiscreteScenarioSet::set_config( Configuration * config )
{
 if( ! config )
  return;  // Nothing to do with null config

 // Pattern 1: SimpleConfiguration<int> - baseline method (top poolSize by
 // weight)
 auto * simple_poolSize = dynamic_cast<SimpleConfiguration<int> *>( config );
 if( simple_poolSize ) {
  ScenarioIndex pool_size = simple_poolSize->f_value;
  if( pool_size > 0 )
   // Baseline method: just set poolSize, no CFL config needed
   poolSize = pool_size;
  return;
  }

 // Pattern 2: SimpleConfiguration<pair<int, Configuration*>> where
 // Configuration* is BlockSolverConfig*
 auto * poolSize_solver_pair =
     dynamic_cast<SimpleConfiguration<std::pair<int, Configuration *>> *>(
         config );
 if( poolSize_solver_pair ) {
  ScenarioIndex pool_size = poolSize_solver_pair->f_value.first;
  Configuration * inner_config = poolSize_solver_pair->f_value.second;

  if( pool_size > 0 && inner_config ) {
   // Check if inner_config is actually a BlockSolverConfig*
   auto * solver_config = dynamic_cast<BlockSolverConfig *>( inner_config );
   if( solver_config ) {
    // Generate simple BlockConfig with poolSize
    auto * block_cfg = generate_default_cfl_config( pool_size );
    // Use advanced scenario reduction with provided solver
    set_config( block_cfg, solver_config, pool_size );
    // Clean up the temporary BlockConfig (set_config clones it)
    delete block_cfg;
   }
  }
  return;
 }

 // Pattern 3: SimpleConfiguration<pair<Configuration*, Configuration*>> where
 // first is BlockConfig*, second is BlockSolverConfig*
 auto * config_pair = dynamic_cast<
     SimpleConfiguration<std::pair<Configuration *, Configuration *>> *>(
     config );
 if( config_pair ) {
  Configuration * first_config = config_pair->f_value.first;
  Configuration * second_config = config_pair->f_value.second;

  if( first_config && second_config ) {
   // Check if first is BlockConfig* and second is BlockSolverConfig*
   auto * block_config = dynamic_cast<BlockConfig *>( first_config );
   auto * solver_config = dynamic_cast<BlockSolverConfig *>( second_config );

   if( block_config && solver_config ) {
    // Extract poolSize from the BlockConfig's extra configuration
    if( !block_config->f_extra_Configuration ) {
     throw std::runtime_error(
         "No extra configuration found containing poolSize parameter" );
    }

    auto * simple_config = dynamic_cast<SimpleConfiguration<int> *>(
        block_config->f_extra_Configuration );
    if( !simple_config ) {
     throw std::runtime_error(
         "Extra configuration is not a SimpleConfiguration<int>" );
    }

    ScenarioIndex pool_size = simple_config->f_value;

    // Validate the extracted poolSize
    validate_poolSize_value( pool_size, nbScenarios );

    // Pass the original BlockConfig (set_config will clone it)
    set_config( block_config, solver_config, pool_size );
   }
  }
  return;
 }


 // If we reach here, the configuration type is not supported
 throw std::invalid_argument(
     "Unsupported configuration type for DiscreteScenarioSet::set_config(). "
     "Supported patterns are:\n"
     "1. SimpleConfiguration<int> - baseline method (top poolSize scenarios by "
     "weight)\n"
     "2. SimpleConfiguration<pair<int, Configuration*>> - advanced with "
     "generated BlockConfig\n"
     "3. SimpleConfiguration<pair<Configuration*, Configuration*>> - full "
     "advanced configuration\n"
     "Note: The ell parameter should be set via set_ell() method or netCDF "
     "deserialization." );
}


/*--------------------------------------------------------------------------*/
/*--------------------- CLASS DiscreteScenarioSet --------------------------*/
/*--------------------------------------------------------------------------*/

void DiscreteScenarioSet::deserialize( const netCDF::NcGroup & group )
{
 // ScenarioGenerator::deserialize is pure virtual, so we start directly here

 // Reset state
 currentScenarioIndex = 0;
 poolSize = 0;
 sumPoolWeights = 0.0;
 is_initialized = false;

 // Clear existing data
 scenarioSet.resize( boost::extents[0][0] );
 setWeights.clear();
 scenarioIndexes.clear();
 universeIndexes.clear();

 // Clean up any existing configuration
 if( f_block_config ) {
  delete f_block_config;
  f_block_config = nullptr;
  }

 // Deserialize mandatory dimensions
 ::deserialize_dim( group, "NumberScenarios", nbScenarios, false );
 ::deserialize_dim( group, "ScenarioSize", scenarioSize, false );

 if( nbScenarios == 0 )
  throw std::invalid_argument( "NumberScenarios must be positive" );

 if( scenarioSize == 0 )
  throw std::invalid_argument( "ScenarioSize must be positive" );

 // Deserialize mandatory scenarios data
 scenarioSet.resize( boost::extents[nbScenarios][scenarioSize] );
 std::vector<std::size_t> sizes = { nbScenarios, scenarioSize };
 ::deserialize( group, "Scenarios", sizes, scenarioSet, true, false );

 // Reserve space for setWeights before deserializing
 setWeights.reserve( nbScenarios );

 // Deserialize probabilities (optional, default to uniform)
 bool probsLoaded =
     ::deserialize( group, "PoolWeights", nbScenarios, setWeights );

 if( ! probsLoaded )
  // No probabilities in file, create uniform weights
  setWeights.assign( nbScenarios, 1.0 / nbScenarios );
 else
  if( setWeights.size() != nbScenarios )
   // Probabilities were loaded but have wrong size - this is an error
   throw std::invalid_argument( "PoolWeights size (" +
				std::to_string( setWeights.size() ) +
				") does not match NumberScenarios (" +
				std::to_string( nbScenarios ) + ")" );

 // Validate probabilities sum to approximately 1.0
 double sum = std::accumulate( setWeights.begin(), setWeights.end(), 0.0 );
 if( std::abs( sum - 1.0 ) > 1e-6 )
  throw std::invalid_argument( "Scenario probabilities sum to "
			       + std::to_string( sum ) + "rather than 1.0" );

 // Deserialize scenario reduction configuration (optional)
 try {
  netCDF::NcGroup cfgGroup = group.getGroup( "ScenarioReductionConfig" );
  if( !cfgGroup.isNull() ) {
   // Try to read poolSize (optional)
   ScenarioIndex pool_size = 0;
   bool has_poolSize_variable = false;
   try {
    auto poolSizeVar = cfgGroup.getVar( "poolSize" );
    if( !poolSizeVar.isNull() ) {
     poolSizeVar.getVar( &pool_size );
     has_poolSize_variable = true;
     // Note: poolSize validation will happen later during init_representative_pool
    }
   }
   catch ( ... ) {
    // poolSize not found, will check BlockConfig for it
   }

   // Try to read ell (optional, default 2.0)
   float ell_value = 2.0f;
   try {
    auto ellVar = cfgGroup.getVar( "ell" );
    if( !ellVar.isNull() ) {
     ellVar.getVar( &ell_value );
    }
   }
   catch ( ... ) {
    // ell not found, use default
   }

   // Store ell as an internal variable
   this->ell = ell_value;

   // Try to deserialize BlockConfig
   BlockConfig * block_cfg = nullptr;
   try {
    auto blockGroup = cfgGroup.getGroup( "BlockConfig" );
    if( !blockGroup.isNull() ) {
     auto * cfg = Configuration::new_Configuration( blockGroup );
     block_cfg = dynamic_cast<BlockConfig *>( cfg );
     if( !block_cfg ) {
      delete cfg;
      throw std::runtime_error( "BlockConfig deserialization failed" );
     }
    }
   }
   catch ( ... ) {
    // No BlockConfig or invalid
    block_cfg = nullptr;
   }

   // Try to deserialize BlockSolverConfig
   BlockSolverConfig * solver_cfg = nullptr;
   try {
    auto solverGroup = cfgGroup.getGroup( "BlockSolverConfig" );
    if( !solverGroup.isNull() ) {
     auto * cfg = Configuration::new_Configuration( solverGroup );
     solver_cfg = dynamic_cast<BlockSolverConfig *>( cfg );
     if( !solver_cfg ) {
      delete cfg;
      throw std::runtime_error( "BlockSolverConfig deserialization failed" );
     }
    }
   }
   catch ( ... ) {
    // No BlockSolverConfig
    solver_cfg = nullptr;
   }

   // Handle BlockConfig and poolSize interaction
   if( block_cfg ) {
    // Check if BlockConfig has extra_Configuration as SimpleConfiguration<int>
    if( block_cfg->f_extra_Configuration ) {
     auto * poolSize_config = dynamic_cast<SimpleConfiguration<int> *>(
         block_cfg->f_extra_Configuration );
     if( poolSize_config ) {
      if( has_poolSize_variable ) {
       // poolSize variable was provided, override the SimpleConfiguration value
       poolSize_config->f_value = pool_size;
      }
      else {
       // No poolSize variable, use SimpleConfiguration<int> value as poolSize
       pool_size = poolSize_config->f_value;
       // Note: poolSize validation will happen later during init_representative_pool
      }
     }
     else if( has_poolSize_variable ) {
      // Extra config exists but is not SimpleConfiguration<int>, and we have poolSize
      delete block_cfg->f_extra_Configuration;
      block_cfg->f_extra_Configuration = new SimpleConfiguration<int>(
          pool_size );
     }
     else {
      // No poolSize variable and extra_config is not SimpleConfiguration<int>
      throw std::runtime_error(
          "poolSize not found: neither as variable nor in BlockConfig's extra_Configuration" );
     }
    }
    else if( has_poolSize_variable ) {
     // No extra_Configuration but poolSize was provided, add it
     block_cfg->f_extra_Configuration = new SimpleConfiguration<int>(
         pool_size );
    }
    else {
     // No poolSize variable and no extra_Configuration
     throw std::runtime_error(
         "poolSize not found: neither as variable nor in BlockConfig" );
    }
   }
   else if( has_poolSize_variable ) {
    // No BlockConfig provided but poolSize was given, create default one with poolSize
    block_cfg = generate_default_cfl_config( pool_size );
    if( !block_cfg ) {
     throw std::runtime_error( "Failed to generate default BlockConfig" );
    }
   }
   else {
    // No BlockConfig and no poolSize variable
    throw std::runtime_error(
        "poolSize not found: must be provided either as variable or in BlockConfig" );
   }

   // Store the configuration
   try {
    set_config( block_cfg, solver_cfg );
    poolSize = pool_size;
   }
   catch ( const std::exception & e ) {
    // Clean up allocated configs on error
    if( block_cfg ) delete block_cfg;
    if( solver_cfg ) delete solver_cfg;
    throw std::runtime_error(
        std::string( "Failed to set scenario reduction config: " ) + e.what() );
   }
  }
 }
 catch ( const std::exception & e ) {
  // Log warning but continue - scenario reduction config is optional
  std::cerr << "Warning: Failed to load ScenarioReductionConfig: " << e.what()
            << std::endl;
 }

 // Lazy canonical init: leave the generator immediately walkable, with the
 // universe set to the full scenario set and the iteration scheduled on the
 // canonical (natural) order. Consumers (e.g., SDDPBlock::deserialize() and
 // SDDPSolver::set_Block() in the default "no pool-size parameter" mode) rely
 // on this so they can iterate the pool right after deserialization without
 // having to call init_*_pool( INFScenario ) explicitly.
 universeIndexes.resize( nbScenarios );
 std::iota( universeIndexes.begin(), universeIndexes.end(), 0 );
 scenarioIndexes = universeIndexes;
 update_pool_weights();
 currentScenarioIndex = 0;
 is_initialized = true;
 }

/*--------------------------------------------------------------------------*/

void DiscreteScenarioSet::serialize( netCDF::NcGroup & group ) const
{
 // Always call base class serialize first
 ScenarioGenerator::serialize( group );

 // Serialize mandatory dimensions
 auto nbScenariosDim = group.addDim( "NumberScenarios", nbScenarios );
 auto scenarioSizeDim = group.addDim( "ScenarioSize", scenarioSize );

 // Serialize mandatory scenarios data
 auto scenariosVar = group.addVar( "Scenarios", netCDF::NcDouble(),
                                   { nbScenariosDim, scenarioSizeDim } );
 scenariosVar.putVar( scenarioSet.data() );

 // Serialize probabilities
 auto probsVar =
     group.addVar( "PoolWeights", netCDF::NcDouble(), nbScenariosDim );
 probsVar.putVar( setWeights.data() );

 // ScenarioReductionConfig is NOT serialized
 // If scenario reduction functionality is needed after deserialization,
 // the configuration must be set separately using set_config()
}


/*--------------------------------------------------------------------------*/
// Load scenario data from a plain-text input stream

void DiscreteScenarioSet::load( std::istream & input )
{
 static const std::string sre( "DiscreteScenarioSet::load: stream read error" );

 // Read N (number of scenarios) and D (scenario size)
 unsigned int N = 0, D = 0;
 input >> N >> D;
 if( input.fail() )
  throw std::invalid_argument( sre + ": expected N D" );
 if( N == 0 )
  throw std::invalid_argument( sre + ": N must be positive" );
 if( D == 0 )
  throw std::invalid_argument( sre + ": D must be positive" );

 // Read probability weights
 std::vector< double > weights( N );
 for( unsigned int i = 0 ; i < N ; ++i ) {
  input >> weights[ i ];
  if( input.fail() )
   throw std::invalid_argument( sre + ": failed reading weight " +
                                std::to_string( i ) );
  }

 // Validate weights sum to 1.0
 double sum = std::accumulate( weights.begin() , weights.end() , 0.0 );
 if( std::abs( sum - 1.0 ) > 1e-6 )
  throw std::invalid_argument( sre + ": weights sum to " +
                               std::to_string( sum ) + ", expected 1.0" );

 // Read scenario data
 std::vector< std::vector< double > > scenarios( N ,
                                                  std::vector< double >( D ) );
 for( unsigned int i = 0 ; i < N ; ++i )
  for( unsigned int d = 0 ; d < D ; ++d ) {
   input >> scenarios[ i ][ d ];
   if( input.fail() )
    throw std::invalid_argument( sre + ": failed reading scenario " +
                                 std::to_string( i ) + " element " +
                                 std::to_string( d ) );
   }

 // Delegate to load_from_memory() which handles validation and state reset
 load_from_memory( scenarios , weights );
 }

/*--------------------------------------------------------------------------*/
// Initialize pool with randomly selected scenarios from the current universe

void DiscreteScenarioSet::init_random_pool( ScenarioIndex pool_size )
{
 // Clean up the iteration state. We do NOT touch universeIndexes here:
 // a previously installed representative filter persists across
 // init_random_pool() calls until explicitly reset by another call to
 // init_representative_pool().
 scenarioIndexes.clear();
 poolWeights.clear();
 sumPoolWeights = 0.0;
 currentScenarioIndex = 0;
 is_initialized = false;

 // Lazily install the default universe (= full set of scenarios) the
 // first time init_random_pool() is invoked without a prior
 // init_representative_pool() call.
 if( universeIndexes.empty() ) {
  universeIndexes.resize( nbScenarios );
  std::iota( universeIndexes.begin() , universeIndexes.end() , 0 );
  }

 const ScenarioIndex universe_size =
  static_cast< ScenarioIndex >( universeIndexes.size() );

 // INFScenario (or any value >= universe size) means "use the whole
 // current universe". For finite values strictly greater than the
 // universe size we throw, since this is almost certainly a caller bug
 // (the "give me everything" intent must be expressed explicitly with
 // INFScenario).
 ScenarioIndex effective_size = pool_size;
 if( pool_size == INFScenario )
  effective_size = universe_size;
 else if( pool_size == 0 )
  throw std::invalid_argument(
   "DiscreteScenarioSet::init_random_pool: pool_size must be positive "
   "(use INFScenario to mean 'all scenarios in the universe')." );
 else if( pool_size > universe_size )
  throw std::invalid_argument(
   "DiscreteScenarioSet::init_random_pool: requested pool_size (" +
   std::to_string( pool_size ) + ") exceeds the current universe size ("
   + std::to_string( universe_size ) + "). Use INFScenario to opt into "
   "'all scenarios in the universe', or call init_representative_pool() "
   "with a larger filter (or INFScenario to reset)." );

 poolSize = effective_size;

 // Sample effective_size scenarios from universeIndexes weighted by
 // setWeights[ universeIndexes[ i ] ]. To keep the existing weighted-
 // sampling helper unchanged, we build a local weight vector over the
 // universe (length universe_size) and then translate the resulting
 // positions back into the global scenarioSet indexing.
 if( effective_size > 0 ) {
  std::vector< double > universe_weights;
  universe_weights.reserve( universe_size );
  for( const auto idx : universeIndexes )
   universe_weights.push_back( setWeights[ idx ] );

  std::vector< ScenarioIndex > positions_in_universe;
  generateWeightedRandomSubset( universe_size , effective_size ,
                                universe_weights ,
                                positions_in_universe , rng );

  scenarioIndexes.reserve( effective_size );
  for( const auto p : positions_in_universe )
   scenarioIndexes.push_back( universeIndexes[ p ] );

  poolWeights.reserve( effective_size );
  update_pool_weights();

  is_initialized = true;
 }
}

/*--------------------------------------------------------------------------*/
// Restrict the universe to a representative subset (filter), and
// initialise the current pool to that universe in canonical order.

void DiscreteScenarioSet::init_representative_pool(
					    ScenarioIndex target_pool_size )
{
 // Clean up the iteration state. universeIndexes is going to be
 // replaced wholesale below.
 scenarioIndexes.clear();
 poolWeights.clear();
 sumPoolWeights = 0.0;
 currentScenarioIndex = 0;
 is_initialized = false;

 // INFScenario (or target_pool_size == nbScenarios) means "no
 // restriction": reset the universe to the full set of scenarios in
 // natural order and walk it canonically. This is also the case the
 // "default arg" semantics relies on. For finite values strictly
 // greater than nbScenarios we throw, mirroring init_random_pool's
 // contract: the "give me everything" intent must be expressed
 // explicitly with INFScenario.
 if( target_pool_size == INFScenario ||
     target_pool_size == nbScenarios ) {

  universeIndexes.resize( nbScenarios );
  std::iota( universeIndexes.begin() , universeIndexes.end() , 0 );

  scenarioIndexes = universeIndexes;
  poolSize = nbScenarios;
  poolWeights.reserve( nbScenarios );
  update_pool_weights();

  is_initialized = ! scenarioIndexes.empty();
  return;
  }

 if( target_pool_size == 0 )
  throw std::invalid_argument(
   "DiscreteScenarioSet::init_representative_pool: target_pool_size "
   "must be positive (use INFScenario to mean 'no restriction')." );

 if( target_pool_size > nbScenarios )
  throw std::invalid_argument(
   "DiscreteScenarioSet::init_representative_pool: requested "
   "target_pool_size (" + std::to_string( target_pool_size ) +
   ") exceeds the available number of scenarios (" +
   std::to_string( nbScenarios ) + "). Use INFScenario to opt into "
   "'all scenarios in the universe'." );

 // Check if we have a BlockSolverConfig (BlockConfig alone is not useful
 // without solver)
 if( f_solver_config ) {
  // Step 2: Create a CapacitatedFacilityLocationBlock for scenario selection
  auto cflBlock = std::make_unique<CapacitatedFacilityLocationBlock>();

  // Set up the scenario selection problem parameters
  ScenarioIndex n_scenarios = nbScenarios;
  ScenarioSize scenario_size = scenarioSize;

  // Create CFL problem data
  auto [capacities, fixed_costs, demands] =
      create_cfl_problem_data( n_scenarios, setWeights );

  // Compute the transport cost matrix with ell-powered distances
  auto transport_costs =
      compute_transport_cost_matrix( n_scenarios, scenario_size, this->ell );

  // Load the CFL problem into the block
  cflBlock->load( n_scenarios,  // Number of facilities
                  n_scenarios,  // Number of customers
                  std::move( capacities ), std::move( fixed_costs ),
                  std::move( demands ), std::move( transport_costs ),
                  true,  // unsplittable: each scenario assigned to exactly
                         // one representative
                  target_pool_size  // Maximum number of facilities to open
  );

  // Step 3.1: Apply BlockConfig if it exists, else generate minimal default
  if( f_block_config ) {
   // Trust that the user provided appropriate BlockConfig for CFL
   f_block_config->apply( cflBlock.get() );
  } else {
   // Generate minimal default config (mostly empty, just target_pool_size in
   // extra_Configuration)
   auto * default_config = generate_default_cfl_config( target_pool_size );
   default_config->apply( cflBlock.get() );
   delete default_config;
  }

  // Step 3.2: Generate abstract variables (to read the solution)
  cflBlock->generate_abstract_variables();

  // Clone and apply BlockSolverConfig to register a solver with the Block
  // clone config to allow reuse
  auto config_copy =
      std::unique_ptr<BlockSolverConfig>( f_solver_config->clone() );
  config_copy->apply( cflBlock.get() );

  // Get the registered solver and solve
  if( cflBlock->get_registered_solvers().empty() ) {
   throw std::runtime_error(
       "No solver registered after BlockSolverConfig::apply" );
  }

  auto * solver = cflBlock->get_registered_solvers().front();
  if( !solver ) {
   throw std::runtime_error( "Failed to get solver from block" );
  }

#ifndef NDEBUG
  std::cout << "DEBUG [init_representative_pool]: Using solver: "
            << solver->classname() << std::endl;
#endif

  // Solve the problem
  int status = solver->compute();
  if( status != Solver::kOK ) {
   throw std::runtime_error( "Solver failed with status: " +
                             std::to_string( status ) );
  }

  // Ensure solver solution is written to CFL Block variables
  solver->get_var_solution();

  // Extract the solution from the CFL block (including assignments)
  std::vector<ScenarioIndex> scenario_assignments;
  extract_scenarios_and_assignments_from_cfl_block(
      cflBlock.get(), n_scenarios, scenarioIndexes, scenario_assignments );

  // Update weights using assignment information
  update_pool_weights_with_assignments( scenario_assignments );
 } else {
  // No solver config - use baseline method (select top target_pool_size by
  // weight) BlockConfig alone is not useful for baseline selection
  apply_baseline_selection( target_pool_size );

  // For baseline selection, use simple renormalization
  update_pool_weights();
 }

 // Update poolSize to reflect actual number of selected scenarios
 poolSize = scenarioIndexes.size();

 // Mirror the freshly-selected representative subset into universeIndexes
 // so that subsequent init_random_pool() calls sample from the restricted
 // universe rather than from the full nbScenarios set. The canonical
 // walk of the universe is exactly the current scenarioIndexes order.
 universeIndexes = scenarioIndexes;

 // Reset current index
 currentScenarioIndex = 0;

 // Mark the pool as initialized if we have scenarios
 if( !scenarioIndexes.empty() ) {
  is_initialized = true;
 }
}

/*--------------------------------------------------------------------------*/

bool DiscreteScenarioSet::next_scenario( void )
{
 // If no pool is initialized or empty, there are no scenarios to move to
 if( !is_initialized || scenarioIndexes.empty() )
  return false;

 if( currentScenarioIndex < scenarioIndexes.size() - 1 ) {
  // Use prefix increment for efficiency
  ++currentScenarioIndex;
  return true;  // Successfully moved to the next scenario
  }

 return false;  // No more scenario in scenarioPool to move to
 }

/*--------------------------------------------------------------------------*/
/// Concrete implementation of ScenarioGenerator

DiscreteScenarioSet::DiscreteScenarioSet( void )
{
 set_seed( 1337 );  // Default seed for reproducibility
 }

/*--------------------------------------------------------------------------*/
/// Destructor

DiscreteScenarioSet::~DiscreteScenarioSet()
{
 // Reset scenario set
 scenarioSet.resize( boost::extents[0][0] );

 // Clear vectors with shrink_to_fit to release memory back to the system
 scenarioIndexes.clear();
 scenarioIndexes.shrink_to_fit();

 universeIndexes.clear();
 universeIndexes.shrink_to_fit();

 setWeights.clear();
 setWeights.shrink_to_fit();

 poolWeights.clear();
 poolWeights.shrink_to_fit();

 // Release the configuration objects following SMS++ ownership patterns

 // BlockConfig: We own our clone, safe to delete
 if( f_block_config ) {
  delete f_block_config;
  f_block_config = nullptr;
  }

 // BlockSolverConfig: unique_ptr handles cleanup automatically
 }

/*--------------------------------------------------------------------------*/
/*-------------------- HELPER METHODS IMPLEMENTATION -----------------------*/
/*--------------------------------------------------------------------------*/

// Compute the transport cost matrix between scenarios
DiscreteScenarioSet::CMatrix
DiscreteScenarioSet::compute_transport_cost_matrix(
    ScenarioIndex n_scenarios, ScenarioSize scenario_size, float ell ) const
{
 CMatrix transport_costs( boost::extents[n_scenarios][n_scenarios] );

 // Map scenarios to Eigen matrix for easier distance calculations
 Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                Eigen::RowMajor>>
     all_scenarios( scenarioSet.data(), n_scenarios, scenario_size );

 // Compute only upper triangle of the symmetric distance matrix
 for( ScenarioIndex i = 0; i < n_scenarios; ++i ) {
  // Diagonal is always zero (distance from scenario to itself)
  transport_costs[i][i] = 0.0;

  // Compute upper triangle only
  for( ScenarioIndex j = i + 1; j < n_scenarios; ++j ) {
   // Calculate distance between scenarios
   Eigen::VectorXd scenario_i = all_scenarios.row( i );
   Eigen::VectorXd scenario_j = all_scenarios.row( j );

   // Compute once and set both symmetric positions
   double distance = compute_scenario_distance( scenario_i, scenario_j, ell );
   transport_costs[i][j] = distance;
   transport_costs[j][i] = distance;  // Symmetric: d(i,j) = d(j,i)
  }
 }

 return transport_costs;
 }

/*--------------------------------------------------------------------------*/
// Compute distance between two scenarios

double DiscreteScenarioSet::compute_scenario_distance(
    const Eigen::VectorXd & scenario1, const Eigen::VectorXd & scenario2,
    float ell ) const
{
 // Transportation cost = ell-power of the euclidean norm
 return std::pow( ( scenario1 - scenario2 ).norm() , ell );
 }


/*--------------------------------------------------------------------------*/
// Apply baseline selection method

void DiscreteScenarioSet::apply_baseline_selection(
						  ScenarioIndex target_size )
{
 // Create pairs of (index, weight) for sorting
 std::vector<std::pair<ScenarioIndex, double>> indexed_weights;
 indexed_weights.reserve( nbScenarios );

 for( ScenarioIndex i = 0; i < nbScenarios; ++i ) {
  double weight = ( i < setWeights.size() ) ? setWeights[i] : 1.0 / nbScenarios;
  indexed_weights.emplace_back( i, weight );
 }

 // Sort by weight in descending order (highest weights first)
 std::sort(
     indexed_weights.begin(), indexed_weights.end(),
     []( const auto & a, const auto & b ) { return a.second > b.second; } );

 // Select the top target_size scenarios
 scenarioIndexes.clear();
 scenarioIndexes.reserve( target_size );

 for( ScenarioIndex i = 0; i < target_size && i < nbScenarios; ++i ) {
  scenarioIndexes.push_back( indexed_weights[i].first );
 }

 // Sort the selected indices for consistency
 std::sort( scenarioIndexes.begin(), scenarioIndexes.end() );

#ifndef NDEBUG
 std::cout << "DEBUG [apply_baseline_selection]: Selected "
           << scenarioIndexes.size()
           << " scenarios using baseline method (top weights)" << std::endl;
#endif
}

/*--------------------------------------------------------------------------*/
// Update pool weights after scenario selection (for random pool)

void DiscreteScenarioSet::update_pool_weights( void )
{
 // Calculate sum of weights of selected scenarios for normalization
 sumPoolWeights = 0.0;
 for( const auto & idx : scenarioIndexes )
  sumPoolWeights += setWeights[idx];

 // Populate pool weights with renormalized values
 poolWeights.clear();
 poolWeights.reserve( scenarioIndexes.size() );

 if( sumPoolWeights > 0.0 ) {
  for( const auto & idx : scenarioIndexes )
   poolWeights.push_back( setWeights[idx] / sumPoolWeights );
  }
 else {
  // Fallback to uniform if all weights are zero
  double uniform_weight = 1.0 / scenarioIndexes.size();
  poolWeights.assign( scenarioIndexes.size(), uniform_weight );
  }
 }

/*--------------------------------------------------------------------------*/
// Update pool weights with assignment information (for representative pool)

void DiscreteScenarioSet::update_pool_weights_with_assignments(
		    const std::vector<ScenarioIndex> & scenario_assignments )
{
 if( scenario_assignments.size() != nbScenarios )
  throw std::invalid_argument( "scenario_assignments size must equal "
			       "nbScenarios" );

 // Create a map from representative index to its position in scenarioIndexes
 std::unordered_map<ScenarioIndex, size_t> repr_to_pos;
 for( size_t pos = 0; pos < scenarioIndexes.size(); ++pos )
  repr_to_pos[scenarioIndexes[pos]] = pos;

 // Initialize aggregated weights for each representative
 std::vector<double> aggregated_weights( scenarioIndexes.size(), 0.0 );

 // Sum weights of all scenarios assigned to each representative
 for( ScenarioIndex i = 0; i < nbScenarios; ++i ) {
  ScenarioIndex assigned_repr = scenario_assignments[i];

  // Find the position of this representative in scenarioIndexes
  auto it = repr_to_pos.find( assigned_repr );
  if( it != repr_to_pos.end() )
   // Add this scenario's weight to its representative
   aggregated_weights[it->second] += setWeights[i];
  else
   // This shouldn't happen if assignments are correct
   throw std::runtime_error( "Scenario " + std::to_string( i ) +
                             " assigned to non-selected representative " +
                             std::to_string( assigned_repr ) );
  }

 // Calculate total weight (should be ~1.0 if all scenarios are assigned, 0 in
 // default uniform case)
 sumPoolWeights = std::accumulate( aggregated_weights.begin() ,
                                   aggregated_weights.end() , 0.0 );

 // Populate pool weights with aggregated and normalized values
 poolWeights.clear();
 poolWeights.reserve( scenarioIndexes.size() );

 if( sumPoolWeights > 0.0 ) {
  for( const auto & weight : aggregated_weights )
   poolWeights.push_back( weight / sumPoolWeights );
  }
 else {
  // Fallback to uniform if all weights are zero
  double uniform_weight = 1.0 / scenarioIndexes.size();
  poolWeights.assign( scenarioIndexes.size(), uniform_weight );
  }

#ifndef NDEBUG
 std::cout
     << "DEBUG [update_pool_weights_with_assignments]: Aggregated weights:"
     << std::endl;
 for( size_t i = 0; i < scenarioIndexes.size(); ++i )
  std::cout << "  Representative " << scenarioIndexes[i]
            << " has aggregated weight " << poolWeights[i] << std::endl;
#endif
 }

/*--------------------------------------------------------------------------*/
/*------------------ End file DiscreteScenarioSet.cpp ----------------------*/
/*--------------------------------------------------------------------------*/