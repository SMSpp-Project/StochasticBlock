/*--------------------------------------------------------------------------*/
/*-------------------- File tests_BendersBFunction.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Main for testing StochasticBlock.
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Rafael Durbano Lobato
 */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <StochasticBlock.h>

#include <algorithm>
#include <random>
#include <vector>

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
  return int_data;
 }

 template< class T >
 typename std::enable_if< std::is_same_v< T , double > ,
                          std::vector< T > & >::type get_data() {
  return dbl_data;
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
 return Block::get_method< typename FunctionType< S , T >::type >
  ( "DummyBlock::set_data" );
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
 return set;
}

template<>
Range build( int size , int total_size ) {
 assert( size <= total_size );
 std::uniform_int_distribution< int > begin_dist( 0 , total_size - size );
 auto begin = begin_dist( random_engine );
 return Range( begin , begin + size );
}

/*--------------------------------------------------------------------------*/

template< class T >
T build_sequential( int size , int offset = 0 );

template<>
Subset build_sequential( int size , int offset ) {
 Subset set( size );
 std::iota( set.begin() , set.end() , offset );
 return set;
}

template<>
Range build_sequential( int size , int offset ) {
 return Range( offset , size + offset );
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
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main() {

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
}
