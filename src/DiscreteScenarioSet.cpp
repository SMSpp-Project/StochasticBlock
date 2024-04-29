/*--------------------------------------------------------------------------*/
/*-------------------- File DiscreteScenarioSet.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Class DiscreteScenarioSet that is an implementation of ScenarioGenerator 
 * suited to the case where the input distribution is contained in a netCDF file
 * as a collection of vectors. Associated with the header file 
 * DiscreteScenarioSet.h
 *
 * \author Antonio Frangioni \n Dipartimento di Informatica \n Universita' di
 *         Pisa \n
 *
 * \author Benoit Tran \n Dipartimento di Informatica \n Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DiscreteScenarioSet.h"
#include "Kmeans.cpp"

#include <stdexcept> // For standard exceptions

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it {
/*--------------------------------------------------------------------------*/
/*--------------------- CLASS DiscreteScenarioSet --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/

/// Concrete implementation of ScenarioGenerator
DiscreteScenarioSet::DiscreteScenarioSet(){
        set_seed(1337);
}

/// Destructor
DiscreteScenarioSet::~DiscreteScenarioSet() {
}
    
/// Implementation for setting the seed of the pseudo-random number generator
void DiscreteScenarioSet::set_seed(unsigned long seed) {
    rng.seed(seed);
}


/// Function to select a random subset from the input scenario pool
/** The function init_random_pool selects randomly a few scenarios among the 
 * ones that were deserialized from the input. It saves a subset of 
 * ScenarioIndex-es into an internal variable called scenarioIndexes.*/
void DiscreteScenarioSet::init_random_pool(ScenarioIndex size) {
    update_poolSize(size);

    // Update in-place ScenarioIndexes
    generateRandomSubset(nbScenarios, size); 
}

/// Function to compute a finite set of representative scenario
/** The function init_representative_pool(ScenarioIndex size) computes a set of 
 * scenarios that approximates the input set of scenarios according to some 
 * optimization criterium. 
 * 
 * By default, a native (and naive) implementation of kmeans clustering is used.
 * 
 * First transform every scenario, contained in a boost::multi_array<double,2>,
 * into an Eigen::VectorXd. Transformed using Eigen::Map, 
 * */
void DiscreteScenarioSet::init_representative_pool(ScenarioIndex size) {
    update_poolSize(size);

    // Convert every input scenario into an Eigen::VectorXd
    PoolMap eigenSet;
    eigenSet.reserve(nbScenarios);
    for (size_t i = 0; i < scenarioSet.shape()[0]; i++)
    {
        eigenSet.emplace_back(Eigen::Map<Eigen::VectorXd>(&scenarioSet[i][0],
         scenarioSize));
    }

    // Updates first rows of scenarioSet with contains the cluster barycenters
    kMeans(size, eigenSet);

    // The representative scenarios are now the first rows of scenarioSet
    scenarioIndexes.resize(size);  
    std::iota(scenarioIndexes.begin(), scenarioIndexes.end(), 0);
}


/// Function for retrieving the current scenario.
/** Checks that the internal variable currentScenarioIndex is within bounds,
 * then converts the currentScenarioIndex-th row of the scenario pool as a 
 * Scenario, that is as a std::span< const double >. */
ScenarioGenerator::Scenario DiscreteScenarioSet::get_current_scenario(void) {
    if (currentScenarioIndex >= poolSize) {
        throw std::out_of_range("Current scenario index is out of range.");
    }
    return Scenario( &scenarioSet[currentScenarioIndex][0],
     get_scenario_size() );
}

double DiscreteScenarioSet::get_current_scenario_probability(void) {
    if (currentScenarioIndex >= poolSize) {
        throw std::out_of_range("Current scenario index is out of range.");
    }
    return scenarioProbabilities[currentScenarioIndex];
}

bool DiscreteScenarioSet::next_scenario(void) {
    if (currentScenarioIndex < poolSize - 1) {
        currentScenarioIndex++;
        return true; // Successfully moved to the next scenario
    }
    return false; // No more scenario in scenarioPool to move to
}

ScenarioGenerator::ScenarioSize DiscreteScenarioSet::get_scenario_size(void) {
    // Implementation for retrieving the size of a scenario
    return scenarioSize;
}

/// deserialize Scenarios from a netCDF::NcGroup
/** Save the Scenarios contained in a netCDF::NcGroup into the internal variable 
 * boost::multi_array< double, 2 > scenarioPool. Optionally, if weights for the 
 * scenarios are provided then save them as well into a std::vector< double >, 
 * otherwise assume uniform weights. */
void DiscreteScenarioSet::deserialize( const netCDF::NcGroup& group ) {
    // Compute the two dimensions of the scenarioPool
    ::SMSpp_di_unipi_it::deserialize_dim( group , "NumberScenarios" ,
                                    nbScenarios , false );

    ::SMSpp_di_unipi_it::deserialize_dim( group , "ScenarioSize" ,
                                        scenarioSize , false );

    // Deserialize the Scenarios inside the scenarioPool
    scenarioSet.resize(boost::extents[nbScenarios][scenarioSize]);
    ::SMSpp_di_unipi_it::deserialize( group, "Scenarios" , scenarioSet, true,
     false);

    // If weights are not present, assume uniform weights
    if ( !::SMSpp_di_unipi_it::deserialize( group , "ScenarioProbabilities" ,
     scenarioSize, scenarioProbabilities )){
        scenarioProbabilities.resize(scenarioSize, 1.0 / scenarioSize);
    }
}

/*--------------------------------------------------------------------------*/
/*-------------------------- FACTORY MANAGEMENT ----------------------------*/
/*--------------------------------------------------------------------------*/
SMSpp_insert_in_factory_cpp_0(DiscreteScenarioSet);

} // namespace SMSpp_di_unipi_it

/*--------------------------------------------------------------------------*/
/*------------------ End file DiscreteScenarioSet.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
