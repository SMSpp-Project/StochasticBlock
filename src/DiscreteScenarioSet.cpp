#include "DiscreteScenarioSet.h"
#include <stdexcept> // For standard exceptions

namespace SMSpp_di_unipi_it {

// Information is contained in a netCDF file, by default will assume that the
// user want to generate a single scenario. Optionally, the user can precise the
// desired sample size with the second variable "size".
DiscreteScenarioSet::DiscreteScenarioSet(const std::string& file){
        ScenarioGenerator::deserialize(file);
        set_seed(1337);
}

DiscreteScenarioSet::~DiscreteScenarioSet() {
    // Destructor implementation 
}

void DiscreteScenarioSet::set_seed(unsigned long seed) {
    // Implementation for setting the seed of the pseudo-random number generator
    this->seed = seed;
}


void DiscreteScenarioSet::init_random_pool(ScenarioIndex size) {
    // Implementation for initializing a scenarioPool
    update_sampleSize(size); // also asserts that size < nbScenarios

    // The pool is made of size Scenario-s among the nbScenarios in scenarioPool
    generateRandomSubset(seed, nbScenarios, size); 
}

void DiscreteScenarioSet::init_representative_pool(ScenarioIndex size) {
    /* Implementation for initializing a representative scenarioPool.

    I wanted to do kmeans clustering. Native implementation in kmeans.cpp (not
    commited) or should I use OpenCV library for optimized implementation?
    */
   update_sampleSize(size); // also asserts that size < NbScenarios
}

ScenarioGenerator::Scenario DiscreteScenarioSet::get_current_scenario(void) {
    /* 
    Implementation for retrieving the current scenario. Ensure that
    currentScenarioIndex is within bounds.
    */
    if (currentScenarioIndex >= sampleSize) {
        throw std::out_of_range("Current scenario index is out of range.");
    }

    // Query the current scenario in the buffer scenarioData
    scenarioPool.getVar({currentScenarioIndex, 0}, {1, scenarioSize}, scenarioData.data());

    // Convert the current scenario to std::span and return
    return std::span< const double >( scenarioData );
}

double DiscreteScenarioSet::get_current_scenario_probability(void) {
    // Implementation for retrieving the probability of the current scenario
    // Ensure that currentScenarioIndex is within bounds
    if (currentScenarioIndex >= scenarioProbabilities.size()) {
        throw std::out_of_range("Current scenario index is out of range.");
    }
    return scenarioProbabilities[currentScenarioIndex];
}

bool DiscreteScenarioSet::next_scenario(void) {
    // Implementation for moving to the next scenario Check if there is a next
    // scenario and update currentScenarioIndex accordingly
    if (currentScenarioIndex + 1 < sampleSize) {
        currentScenarioIndex++;
        return true; // Successfully moved to the next scenario
    }
    return false; // No more scenario in scenarioPool to move to
}

ScenarioGenerator::ScenarioSize DiscreteScenarioSet::get_scenario_size(void) {
    // Implementation for retrieving the size of a scenario
    return scenarioSize;
}

void DiscreteScenarioSet::deserialize(const netCDF::NcGroup& group) {
    // Implementation for deserializing a DiscreteScenarioSet from a
    // netCDF::NcGroup 

    ::SMSpp_di_unipi_it::deserialize_dim( group , "NumberScenarios" ,
                                    nbScenarios , false );

    ::SMSpp_di_unipi_it::deserialize_dim( group , "ScenarioSize" ,
                                        scenarioSize , false );

    /*
        NEED DESIGN DISCUSSION WITH ANTONIO, (see todo.md not comitted)
        deserialize_scenarios( group );
    */
   scenarioPool = group.getVar("Scenarios");


    // If weights are not present, assume uniform weights
    try{
    ::SMSpp_di_unipi_it::deserialize_dim( group , "ScenarioProbabilities" , scenarioProbabilities, true );
    } catch (const std::invalid_argument &e){
        scenarioProbabilities.resize(scenarioSize, 1.0 / scenarioSize);
    }
}

/* 
    NEED DESIGN DISCUSSION WITH ANTONIO, (see todo.md not comitted)

    // Internal method implementations

    void DiscreteScenarioSet::deserialize_scenarios( const netCDF::NcGroup & group ) {
        // Gather the input scenarios in the scenarioPool
        scenarioPool.resize( nbScenarios , std::vector< double >( scenarioSize ) );

        auto scenarios_var = group.getVar( "Scenarios" );

        if( scenarios_var.isNull() )
        throw( std::invalid_argument
                ( "ScenarioSet::deserialize_scenarios: 'Scenarios' "
                    "variable has not been provided." ) );

        auto dims = scenarios_var.getDims();

        if( ( dims.size() != 2 ) || ( dims[ 0 ].getSize() != nbScenarios ) ||
            ( dims[ 1 ].getSize() != scenarioSize ) )

        throw( std::logic_error
                ( "ScenarioSet::deserialize_scenarios: 'Scenarios' must be a two-"
                    "dimensional array whose first and second dimensions have sizes "
                    "'NumberScenarios' and 'ScenarioSize', respectively." ) );
        for( decltype(scenarioPool)::size_type i = 0 ; i < scenarioPool.size() ; ++i )
        scenarios_var.getVar( { i , 0 } , { 1 , scenarioPool[ i ].size() } ,
                                scenarioPool[ i ].data() );
    }
*/

// Additional functions


// Factory registration. Why doesn't it work?
// SMSpp_insert_in_factory_cpp_0("DiscreteScenarioSet");

} // namespace SMSpp_di_unipi_it
