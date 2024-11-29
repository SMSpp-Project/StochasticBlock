/*--------------------------------------------------------------------------*/
/*------------------- File ScenarioReductorBlock.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the ScenarioReductorBlock class.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 * 
 * \author Romain Pujol \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 * 
 * \author Benoît Tran \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Romain Pujol and Benoît Tran
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ScenarioReductorBlock.h"

#include <iomanip>

/*--------------------------------------------------------------------------*/
/*--------------------------------- MACROS ---------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 #define CHECK_DS 0
 /* Perform long and costly checks on the data structures representing the
  * abstract and the physical representations agree. */
#else
 #define CHECK_DS 0
 // never change this
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- TYPES -----------------------------------*/
/*--------------------------------------------------------------------------*/

using Index = Block::Index;
using c_Index = Block::c_Index;

using Range = Block::Range;
using c_Range = Block::c_Range;

using Subset = Block::Subset;
using c_Subset = Block::c_Subset;

using TNumber = ScenarioReductorBlock::TNumber;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

// static constexpr auto dNAN = std::numeric_limits< double >::quiet_NaN();

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register ScenarioReductorBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( ScenarioReductorBlock );

// register ScenarioReductorSolution to the Solution factory

// SMSpp_insert_in_factory_cpp_0( ScenarioReductorSolution );

/*--------------------------------------------------------------------------*/
/*------------------- METHODS OF ScenarioReductorBlock ---------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------ OTHER INITIALIZATIONS ---------------------------*/

// load from memory by copying
void ScenarioReductorBlock::load( Index n,
                                  Index m,
                                  Index dim,
                                  c_Mat_ANumber & atoms,
                                  c_Vec_WNumber & weights,
                                  Index k )
{
 // Sanity checks 

 if( atoms.size() != n )
  throw( std::invalid_argument( "Vector of atoms of the wrong size" ) );

 if( weights.size() != n )
  throw( std::invalid_argument( "Vector of weights of the wrong size" ) );

 for( int i{0} ; i < n ; i++ )
 {
  if( atoms[ i ].size() != dim )
    throw( std::invalid_argument( "The given dimension is wrong" ) );
 }
 
 // call load( , , && , && , ) on newly constructed copies
 load( n,
       m,
       dim,
       Mat_ANumber( atoms ),
       Vec_WNumber( weights ),
       k );

 } // end( ScenarioReductorBlock::load( memory ) )

/*--------------------------------------------------------------------------*/
// load from memory by moving
 void ScenarioReductorBlock::load(  Index n,
                                    Index m,
                                    Index dim,
                                    c_Mat_ANumber && atoms,
                                    c_Vec_WNumber && weights,
                                    Index k )
{
// Sanity checks 

 if( atoms.size() != n )
  throw( std::invalid_argument( "Vector of atoms of the wrong size" ) );

 if( weights.size() != n )
  throw( std::invalid_argument( "Vector of weights of the wrong size" ) );

 for( int i{0} ; i < n ; i++ )
 {
  if( atoms[ i ].size() != dim )
    throw( std::invalid_argument( "The given dimension is wrong" ) );
 }

 // Erase previous instance, if any
 if( get_N() ) {
  // guts_of_destructor(); TODO fixme
 }

 // Initialize physical fields of ScenarioReductorBlock
 this->N = n;
 this->M = m;
 this->dim = dim;
 this->k = k;
 
 atomsP = std::move( atoms );
 weightsP = std::move( weightsP );
 
 generate_abstract_variables();

 // Declare modification
 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );
}  // end( ScenarioReductorBlock::load( memory ) )

/*--------------------------------------------------------------------------*/
void ScenarioReductorBlock::load( std::istream & input , char frmt )
{
 // erase previous instance, if any
 if( get_N() ) {
  // guts_of_destructor(); TODO fixme
 }

 // read problem data
 Index n;
 if( ! ( input >> eatcomments >> n ) )
  throw( std::invalid_argument( "error reading number of atoms of P" ) );

 this->N = n;
 // atomsP.resize( n ); TODO fixme
 weightsP.resize( n );

 // read the desired output size
 Index m;
 if( ! ( input >> eatcomments >> m ) )
  throw( std::invalid_argument( "error reading number of desired atoms of Q" ) );

 this->M = m;

 // read dimension
 Index D;
 if( ! ( input >> eatcomments >> D ) )
  throw( std::invalid_argument( "error reading the dimension" ) );
 this->dim = D;
 atomsP.resize( boost::extents[ n ][ dim ] );

 // read Wasserstein order
 Index k;
 if( ! ( input >> eatcomments >> k ) )
  throw( std::invalid_argument( "error reading the Wasserstein order" ) );
 this->k = k;

 // read atoms
 for( Index i = 0 ; i < get_N() ; ++i )
  {
    for( Index d = 0 ; d < get_dim() ; ++d )
    {
      if( ! ( input >> eatcomments >> atomsP[ i ][ d ] ) )
        throw( std::invalid_argument( "error reading the atomsP" ) );
    }
  }

 // read weights
 for( Index i = 0 ; i < get_N() ; ++i )
  if( ! ( input >> eatcomments >> weightsP[ i ] ) )
   throw( std::invalid_argument( "error reading Weights" ) );

 generate_abstract_variables();

 // Modification- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );

 } // end( BinaryKnapsackBlock::load( istream ) )

/*--------------------------------------------------------------------------*/
void ScenarioReductorBlock::deserialize( const netCDF::NcGroup & group )
{
 // erase previous instance, if any
 if( get_N() ) {
  // guts_of_destructor();  TODO fixme
 }

 // read N
 ::deserialize_dim( group , "N" , N );

 // read M
 ::deserialize_dim( group , "M" , M );

 // read dim
 ::deserialize_dim( group , "dim" , dim );

 // read AtomsP
 atomsP.resize( boost::extents[ get_N() ][ get_dim() ] );
 ::deserialize( group , "atomsP" , atomsP , true , false );

 // read weightsP, if not present, assume uniform weights
 if( ! ::deserialize( group , "weightsP" , get_N() , weightsP ) ) {
  weightsP.resize( get_N() , 1.0 / get_N() );
 }

 generate_abstract_variables();

 // call the method of Block
 // inside this the NBModification, the "nuclear option", is issued
 Block::deserialize( group );

 }  // end( ScenarioReductorBlock::deserialize )

/*--------------------------------------------------------------------------*/
/*------------------------- Methods for R3 Blocks --------------------------*/
/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---- METHODS FOR LOADING, PRINTING & SAVING THE ScenarioReductorBlock ----*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*---------------- End File ScenarioReductorBlock.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
