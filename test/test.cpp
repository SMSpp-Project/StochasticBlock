/*--------------------------------------------------------------------------*/
/*---------------------------- File tests.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Main for testing StochasticBlock.
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 * 
 * \author Benoît Tran \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Rafael Durbano Lobato, Benoît Tran
 */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <StochasticBlock.h>
#include <DiscreteScenarioSet.h>  // For concrete testing

#include <random>
#include <iostream>
#include <cstdio>  // For std::remove

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using Subset = Block::Subset;
using Range = Block::Range;

/*--------------------------------------------------------------------------*/
/*-------------------------- GLOBAL VARIABLES ------------------------------*/
/*--------------------------------------------------------------------------*/

std::mt19937 random_engine;

/*--------------------------------------------------------------------------*/
/*--------------------------- AUXILIARY TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

class DummyBlock : public Block {

public:

 DummyBlock( Block * f_block = nullptr ) : Block( f_block ) {}

 DummyBlock( std::size_t int_size , std::size_t dbl_size ) {
  int_data.resize( int_size );
  std::iota( int_data.begin() , int_data.end() , 0 );

  dbl_data.resize( dbl_size );
  std::iota( dbl_data.begin() , dbl_data.end() , 0.0 );
 }

 template< class T >
 typename std::enable_if< std::is_same_v< T , int > ,
                          std::vector< T > & >::type get_data() {
  return( int_data );
 }

 template< class T >
 typename std::enable_if< std::is_same_v< T , double > ,
                          std::vector< T > & >::type get_data() {
  return( dbl_data );
 }

 template< class T >
 void set_data( typename std::vector< T >::const_iterator values ,
                Subset && subset ,
                bool ordered = false ,
                c_ModParam issuePMod = eNoBlck ,
                c_ModParam issueAMod = eNoBlck ) {
  auto & data = get_data< T >();
  for( auto i : subset ) {
   data[ i ] = *values;
   std::advance( values , 1 );
  }
 }

 template< class T >
 void set_data( typename std::vector< T >::const_iterator values ,
                Range rng = Range( 0, Inf< Index >() ) ,
                c_ModParam issuePMod = eNoBlck ,
                c_ModParam issueAMod = eNoBlck ) {
  auto & data = get_data< T >();
  decltype( rng.second ) size = data.size();
  rng.second = std::min( size , rng.second );
  for( std::size_t i = rng.first ; i < rng.second ; ++i ) {
   data[ i ] = *values;
   std::advance( values , 1 );
  }
 }

 static void static_initialization() {
  register_method< DummyBlock , MF_int_it , Subset && , const bool >
   ( "DummyBlock::set_data" , & DummyBlock::set_data< int > );

  register_method< DummyBlock , MF_int_it , Range >
   ( "DummyBlock::set_data" , & DummyBlock::set_data< int > );

  register_method< DummyBlock , MF_dbl_it , Subset && , const bool >
   ( "DummyBlock::set_data" , & DummyBlock::set_data< double > );

  register_method< DummyBlock , MF_dbl_it , Range >
   ( "DummyBlock::set_data" , & DummyBlock::set_data< double > );
 }

protected:

 void load( std::istream & input , char frmt ) override {}

private:
 std::vector< int > int_data;
 std::vector< double > dbl_data;
 SMSpp_insert_in_factory_h;
};

SMSpp_insert_in_factory_cpp_1( DummyBlock );

/*--------------------------------------------------------------------------*/
/*------------------------- AUXILIARY FUNCTIONS ----------------------------*/
/*--------------------------------------------------------------------------*/

template< class T >
struct Iter;

template<>
struct Iter< int > {
 using type = Block::MF_int_it;
};

template<>
struct Iter< double > {
 using type = Block::MF_dbl_it;
};

template< class S , class T >
struct FunctionType;

template< class T >
struct FunctionType< Subset , T > {
 using type = Block::FunctionType< typename Iter< T >::type , Subset && , bool >;
};

template< class T >
struct FunctionType< Range , T > {
 using type = Block::FunctionType< typename Iter< T >::type , Range >;
};

template< class S , class T >
auto get_method() {
 return( Block::get_method< typename FunctionType< S , T >::type >
  ( "DummyBlock::set_data" ) );
}

/*--------------------------------------------------------------------------*/

template< class T >
T build( int size , int total_size );

template<>
Subset build( int size , int total_size ) {
 assert( size <= total_size );
 Subset set( total_size );
 std::iota( set.begin() , set.end() , 0 );
 std::shuffle( set.begin() , set.end() , random_engine );
 set.resize( size );
 std::sort( set.begin() , set.end() );
 return( set );
}

template<>
Range build( int size , int total_size ) {
 assert( size <= total_size );
 std::uniform_int_distribution< int > begin_dist( 0 , total_size - size );
 auto begin = begin_dist( random_engine );
 return( Range( begin , begin + size ) );
}

/*--------------------------------------------------------------------------*/

template< class T >
T build_sequential( int size , int offset = 0 );

template<>
Subset build_sequential( int size , int offset ) {
 Subset set( size );
 std::iota( set.begin() , set.end() , offset );
 return( set );
}

template<>
Range build_sequential( int size , int offset ) {
 return( Range( offset , size + offset ) );
}

/*--------------------------------------------------------------------------*/

template< class SetTo , class T >
void check( SetTo set_to , const std::vector< double > & data ,
            const std::vector< T > & block_data , double offset );

// we assume set_to is ordered
template< class T >
void check( Subset set_to , const std::vector< double > & data ,
            const std::vector< T > & block_data ) {
 std::size_t j = 0;
 for( std::size_t i = 0 ; i < block_data.size() ; ++i ) {
  if( j < set_to.size() && set_to[ j ] == i ) {
   assert( std::size_t( block_data[ set_to[ j ] ] ) == data[ j ] );
   ++j;
  }
  else
   assert( std::size_t( block_data[ i ] ) == i  );
 }
 assert( j == set_to.size() );
}

template< class T >
void check( Range set_to , const std::vector< double > & data ,
            const std::vector< T > & block_data ) {

 for( std::size_t i = 0 ; i < set_to.first ; ++i ) {
  assert( std::size_t( block_data[ i ] ) == i );
 }

 for( std::size_t i = set_to.first ; i < set_to.second ; ++i ) {
  assert( std::size_t( block_data[ i ] ) == data[ i - set_to.first ] );
 }

 for( std::size_t i = set_to.second ; i < data.size() ; ++i ) {
  assert( std::size_t( block_data[ i ] ) == i );
 }
}

/*--------------------------------------------------------------------------*/

// Helper function to create a test netCDF file with scenarios
void create_test_scenario_file(const std::string& filename,
                               int num_scenarios,
                               int scenario_size) {
 // Create a new netCDF file
 netCDF::NcFile dataFile(filename, netCDF::NcFile::replace);
 
 // Add dimensions
 auto nbScenariosDim = dataFile.addDim("NumberScenarios", num_scenarios);
 auto scenarioSizeDim = dataFile.addDim("ScenarioSize", scenario_size);
 
 // Create and fill scenario data
 std::vector<double> scenarios(num_scenarios * scenario_size);
 for(int s = 0; s < num_scenarios; ++s) {
  for(int i = 0; i < scenario_size; ++i) {
   scenarios[s * scenario_size + i] = 100.0 * (s + 1) + i;
  }
 }
 
 // Add Scenarios variable
 auto scenariosVar = dataFile.addVar("Scenarios", netCDF::ncDouble, 
                                     {nbScenariosDim, scenarioSizeDim});
 scenariosVar.putVar(scenarios.data());
 
 // Add weights (uniform for simplicity)
 std::vector<double> weights(num_scenarios, 1.0 / num_scenarios);
 auto weightsVar = dataFile.addVar("poolWeights", netCDF::ncDouble, nbScenariosDim);
 weightsVar.putVar(weights.data());
 
}

/*--------------------------------------------------------------------------*/
/*-------------------------------- TESTS -----------------------------------*/
/*--------------------------------------------------------------------------*/

template< class SetFrom , class SetTo >
void test( std::size_t int_size , std::size_t dbl_size ) {

 auto inner_block = new DummyBlock( int_size , dbl_size );
 StochasticBlock stochastic_block( nullptr , inner_block );

 std::uniform_int_distribution< int > uniform_dist_int( 0 , int_size );
 std::uniform_int_distribution< int > uniform_dist_dbl( 0 , dbl_size );

 std::size_t scenario_int_size = uniform_dist_int( random_engine );
 std::size_t scenario_dbl_size = uniform_dist_dbl( random_engine );

 SetFrom set_from_int = build_sequential< SetFrom >( scenario_int_size );
 SetTo set_to_int = build< SetTo >( scenario_int_size , int_size );

 stochastic_block.add_data_mapping
  ( std::make_unique< SimpleDataMapping< SetFrom , SetTo , int > >
    ( get_method< SetTo , int >() , inner_block ,
      set_from_int , set_to_int ) );

 SetFrom set_from_dbl = build_sequential< SetFrom >
  ( scenario_dbl_size , scenario_int_size );
 SetTo set_to_dbl = build< SetTo >( scenario_dbl_size , dbl_size );

 stochastic_block.add_data_mapping
  ( std::make_unique< SimpleDataMapping< SetFrom , SetTo , double > >
    ( get_method< SetTo , double >() , inner_block ,
      set_from_dbl , set_to_dbl ) );

 std::vector< double > int_data( scenario_int_size );
 std::vector< double > dbl_data( scenario_dbl_size );

 for( std::size_t i = 0 ; i < int_data.size() ; ++i )
  int_data[ i ] = 1.0e6 + i;
 for( std::size_t i = 0 ; i < dbl_data.size() ; ++i )
  dbl_data[ i ] = 2.0e6 + i;

 std::vector< double > data( int_data );
 data.insert( data.end() , dbl_data.begin() , dbl_data.end() );

 stochastic_block.set_data( data.begin() );

 auto block_int_data = inner_block->get_data< int >();
 auto block_dbl_data = inner_block->get_data< double >();

 // check

 check( set_to_int , int_data , block_int_data );
 check( set_to_dbl , dbl_data , block_dbl_data );
}

/*--------------------------------------------------------------------------*/

// Control verbosity of test output
bool VERBOSE_TESTS = false;  // Can be set via command-line argument

void test_scenario_generator() {
 if (VERBOSE_TESTS) std::cout << "Testing ScenarioGenerator integration..." << std::endl;
 
 // Part 1: Test with netCDF file and proper DiscreteScenarioSet initialization
 if (VERBOSE_TESTS) std::cout << "  Part 1: Testing with netCDF file..." << std::endl;
 
 // Create a test netCDF file
 const std::string test_file = "test_scenarios.nc";
 const int num_scenarios = 5;
 const int scenario_size = 8;
 create_test_scenario_file(test_file, num_scenarios, scenario_size);
 
 // Create a DiscreteScenarioSet and load from file
 auto* dss = new DiscreteScenarioSet();
 netCDF::NcFile dataFile(test_file, netCDF::NcFile::read);
 dss->deserialize(dataFile);
 
 // Create a DummyBlock and StochasticBlock
 auto inner_block = new DummyBlock(5, 10);
 StochasticBlock stochastic_block(nullptr, inner_block);
 
 // Set up data mappings
 Subset int_from = {0, 1, 2};
 Subset int_to = {0, 2, 4};
 stochastic_block.add_data_mapping(
  std::make_unique<SimpleDataMapping<Subset, Subset, int>>(
   get_method<Subset, int>(), inner_block, int_from, int_to));
 
 Range dbl_from(3, 8);
 Subset dbl_to = {1, 3, 5, 7, 9};
 stochastic_block.add_data_mapping(
  std::make_unique<SimpleDataMapping<Range, Subset, double>>(
   get_method<Subset, double>(), inner_block, dbl_from, dbl_to));
 
 // Test setter/getter for ScenarioGenerator
 stochastic_block.set_scenario_generator(dss);
 auto* retrieved_sg = stochastic_block.get_scenario_generator();
 assert(retrieved_sg != nullptr);
 assert(dynamic_cast<DiscreteScenarioSet*>(retrieved_sg) != nullptr);
 if (VERBOSE_TESTS) std::cout << "    ScenarioGenerator setter/getter work correctly!" << std::endl;
 
 // Initialize a random pool and iterate through scenarios
 dss->init_random_pool(3);  // Select 3 random scenarios
 
 int scenario_count = 0;
 do {
  auto scenario = dss->get_current_scenario();
  double prob = dss->get_current_scenario_probability();
  
  if (VERBOSE_TESTS) {
   std::cout << "    Scenario " << scenario_count << " (prob=" << prob << "): ";
   for(size_t i = 0; i < std::min(size_t(4), scenario.size()); ++i) {
    std::cout << scenario[i] << " ";
   }
   std::cout << "..." << std::endl;
  }
  
  // Apply scenario to StochasticBlock
  stochastic_block.set_data(scenario);
  
  // Verify application
  auto& int_data = inner_block->get_data<int>();
  assert(int_data[0] == scenario[0]);
  assert(int_data[2] == scenario[1]);
  assert(int_data[4] == scenario[2]);
  
  auto& dbl_data = inner_block->get_data<double>();
  assert(dbl_data[1] == scenario[3]);
  assert(dbl_data[3] == scenario[4]);
  
  scenario_count++;
 } while(dss->next_scenario());
 
 assert(scenario_count == 3);  // We selected 3 random scenarios
 if (VERBOSE_TESTS) std::cout << "    Processed " << scenario_count << " scenarios successfully!" << std::endl;
 
 // Clean up test file
 std::remove(test_file.c_str());
 
 // Part 2: Test direct set_data and new set_data(Scenario) overload
 if (VERBOSE_TESTS) std::cout << "  Part 2: Testing direct set_data methods..." << std::endl;
 
 // Reuse the same inner_block and stochastic_block from Part 1
 // Test with manually created scenario data
 const int manual_scenarios = 3;
 for( int s = 0; s < manual_scenarios; ++s ) {
  std::vector<double> scenario_data(scenario_size);
  for( int i = 0; i < scenario_size; ++i ) {
   scenario_data[i] = 100.0 * (s + 1) + i;  // Scenario 1: 100-107, Scenario 2: 200-207, etc.
  }
  
  if (VERBOSE_TESTS && s == 0) {
   std::cout << "    Testing vector set_data with scenario: ";
   for(size_t i = 0; i < std::min(size_t(4), scenario_data.size()); ++i) {
    std::cout << scenario_data[i] << " ";
   }
   std::cout << "..." << std::endl;
  }
  
  // Apply scenario using set_data
  stochastic_block.set_data( scenario_data );
  
  // Verify the data was applied correctly
  auto& int_data_check = inner_block->get_data<int>();
  auto& dbl_data_check = inner_block->get_data<double>();
  
  // Check int mappings (positions 0,2,4 should have scenario values 0,1,2)
  assert( int_data_check[0] == scenario_data[0] );
  assert( int_data_check[2] == scenario_data[1] );
  assert( int_data_check[4] == scenario_data[2] );
  
  // Check double mappings (positions 1,3,5,7,9 should have scenario values 3-7)
  assert( dbl_data_check[1] == scenario_data[3] );
  assert( dbl_data_check[3] == scenario_data[4] );
  assert( dbl_data_check[5] == scenario_data[5] );
  assert( dbl_data_check[7] == scenario_data[6] );
  assert( dbl_data_check[9] == scenario_data[7] );
 }
 
 // Test the new set_data(Scenario) overload
 std::vector<double> test_scenario = {10, 20, 30, 40, 50, 60, 70, 80};
 ScenarioGenerator::Scenario scenario_span(test_scenario.data(), test_scenario.size());
 stochastic_block.set_data(scenario_span);
 
 auto& final_int_data = inner_block->get_data<int>();
 auto& final_dbl_data = inner_block->get_data<double>();
 assert( final_int_data[0] == 10 );
 assert( final_int_data[2] == 20 );
 assert( final_int_data[4] == 30 );
 assert( final_dbl_data[1] == 40 );
 assert( final_dbl_data[3] == 50 );
 if (VERBOSE_TESTS) std::cout << "    set_data(Scenario) overload works correctly!" << std::endl;
 
 if (VERBOSE_TESTS) std::cout << "ScenarioGenerator integration test completed successfully!" << std::endl << std::endl;
}

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main(int argc, char* argv[]) {

 // Parse command-line arguments
 for (int i = 1; i < argc; ++i) {
  std::string arg = argv[i];
  if (arg == "--verbose" || arg == "-v") {
   VERBOSE_TESTS = true;
   std::cout << "Verbose mode enabled" << std::endl;
  } else if (arg == "--help" || arg == "-h") {
   std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
   std::cout << "Options:" << std::endl;
   std::cout << "  -v, --verbose    Enable verbose output" << std::endl;
   std::cout << "  -h, --help       Show this help message" << std::endl;
   return 0;
  } else {
   std::cerr << "Unknown argument: " << arg << std::endl;
   std::cerr << "Use --help for usage information" << std::endl;
   return 1;
  }
 }

 std::uniform_int_distribution< int > size_dist( 0 , 20 );

 for( int i = 0 ; i < 10000 ; ++i ) {
  test< Subset , Subset >( size_dist( random_engine ) ,
                           size_dist( random_engine ) );

  test< Subset , Range >( size_dist( random_engine ) ,
                          size_dist( random_engine ) );

  test< Range , Subset >( size_dist( random_engine ) ,
                          size_dist( random_engine ) );

  test< Range , Range >( size_dist( random_engine ) ,
                         size_dist( random_engine ) );
 }
 
 // Unit test for ScenarioGenerator integration
 test_scenario_generator(); 
 
 if (VERBOSE_TESTS) {
  std::cout << "All tests passed successfully!" << std::endl;
 }
}

/*--------------------------------------------------------------------------*/
/*------------------------- End File test.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
