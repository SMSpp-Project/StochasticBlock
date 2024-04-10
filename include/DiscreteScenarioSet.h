/*--------------------------------------------------------------------------*/
/*-------------------- File DiscreteScenarioSet.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 *
 * \author Benoit Tran \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Benoit Tran
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __DiscreteScenarioSet
#define __DiscreteScenarioSet
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ScenarioGenerator.h"

#include <vector>
#include <random>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup
 *  @{ */

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS DiscreteScenarioSet --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
///
/**
 */

class DiscreteScenarioSet : public ScenarioGenerator
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*------------- CONSTRUCTING AND DESTRUCTING StochasticBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing DiscreteScenarioSet
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// constructor
 // Information should be contained in a netCDF file
 DiscreteScenarioSet( const std::string & file );

 virtual ~DiscreteScenarioSet();

 // Implementing pure virtual methods from ScenarioGenerator
 virtual void set_seed( unsigned long seed ) override;

 virtual void init_random_pool( ScenarioIndex sampleSize ) override;

 virtual void init_representative_pool( ScenarioIndex sampleSize ) override;

 virtual Scenario get_current_scenario( void ) override;

 virtual double get_current_scenario_probability( void ) override;

 virtual bool next_scenario( void ) override;

 virtual ScenarioSize get_scenario_size( void ) override;

 virtual void deserialize( const netCDF::NcGroup & group ) override;

 // Additional methods and attributes specific to DiscreteScenarioSet
 // ...

 protected:
 // Internal methods

 /*
     NEED DESIGN DISCUSSION WITH ANTONIO, (see todo.md not committed)
     virtual void deserialize_scenarios( const netCDF::NcGroup & group ) = 0;
 */
 private:
 // Private attributes
 // std::vector<std::vector<double>> scenarioPool; // Container for scenarios
 netCDF::NcVar scenarioPool; // Container for Scenario-s
 ScenarioSize scenarioSize; // Size of a scenario
 std::vector< ScenarioIndex > scenarioIndexes; // Indexes of the pool
 ScenarioIndex currentScenarioIndex = 0; // Current index in the pool
 std::vector< double > scenarioData; // Buffer for storing a scenario
 std::vector< double > scenarioProbabilities; // Probabilities of scenarios
 unsigned int nbScenarios; // Number of different scenarios
 ScenarioIndex sampleSize; // Desired sample size
 unsigned long seed; // Seed for the pseudo-random number generator



 // Helper methods

 void update_sampleSize( ScenarioIndex size )
 {
  if( size > nbScenarios ) {
   throw std::out_of_range(
    "The desired sample size is greater than the number of available number of different scenarios" );
  }
  sampleSize = size;
 }

 // Draw k elements among n using our seed for reproducibility
 std::vector< int >
 generateRandomSubset( unsigned long seed , size_t n , size_t k )
 {
  std::vector< int > elements( n );
  std::vector< int > subset;

  // Fill the vector with a sequence of n elements (e.g., 0 to n-1)
  for( size_t i = 0 ; i < n ; ++i ) {
   elements[ i ] = i;
  }

  // Initialize a random number generator with the given seed
  std::mt19937 rng( seed );

  // Shuffle the elements randomly
  std::shuffle( elements.begin() , elements.end() , rng );

  // Select the first k elements as the subset
  for( size_t i = 0 ; i < k ; ++i ) {
   subset.push_back( elements[ i ] );
  }

  return subset;
 }

 // Macro for the factory
 SMSpp_insert_in_factory_h;
};

} // namespace SMSpp_di_unipi_it

#endif // __DiscreteScenarioSet
