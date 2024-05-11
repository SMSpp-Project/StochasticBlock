/*--------------------------------------------------------------------------*/
/*------------------------- File test_discrete.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Main for testing DiscreteScenarioSet.
 * 
 * Put the compiling option -verbose for additional comments.
 *
 * \author Benoît Tran \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Benoît Tran
 */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/
#include "DiscreteScenarioSet.h"

#include <iostream>
#include <netcdf>
#include <sstream> // to conditionnaly kill the cout

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/


using namespace netCDF;
using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*--------------------------- CoutSuppressor -------------------------------*/
/*--------------------------------------------------------------------------*/


class CoutSuppressor {
public:
    CoutSuppressor(bool verbose) : verbose(verbose) {
        if (!verbose) {
            // Redirect std::cout to a null stream if not verbose
            old_buf = std::cout.rdbuf(null_buf.rdbuf());
        }
    }

    ~CoutSuppressor() {
        if (!verbose) {
            // Restore std::cout's original stream buffer
            std::cout.rdbuf(old_buf);
        }
    }

    template<typename T>
    CoutSuppressor& operator<<(const T& data) {
        if (verbose) {
            std::cout << data;
        }
        return *this;
    }

private:
    bool verbose;
    std::streambuf* old_buf;
    std::stringstream null_buf;
};

/*--------------------------------------------------------------------------*/
/*------------------------- AUXILIARY FUNCTIONS ----------------------------*/
/*--------------------------------------------------------------------------*/

// Function to create a simple netCDF file
void simpleNetCDF(const std::string& filename, bool probas) {
    NcFile dataFile(filename, NcFile::replace);
    NcDim scenarioDim = dataFile.addDim("NumberScenarios", 5);
    NcDim dimensionDim = dataFile.addDim("ScenarioSize", 10);
    NcVar dataVar = dataFile.addVar("Scenarios", ncDouble, {scenarioDim,
     dimensionDim});
    std::vector<double> scenarios(scenarioDim.getSize() * 
     dimensionDim.getSize(), 0.5);

    /* if probVar not updated, it should be of size 0 and then the deserialize 
    of that var should be false then the default should still make it the 
    constant vector 1/5 */
    NcVar probVar = dataFile.addVar("ScenarioProbabilities", ncDouble,
      scenarioDim);
    if (probas)
    {
      std::vector<double> probabilities(scenarioDim.getSize(),
      1.0 / scenarioDim.getSize());
      probVar.putVar(probabilities.data());
    }
    dataVar.putVar(scenarios.data());
}

/*--------------------------------------------------------------------------*/
/*---------------------------- UNIT TESTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

// Unit test of the deserialization
void test_deserialize()
{
    std::cout << "Deserialization test" << std::endl;
    std::string filename = "temp_scenarioData.nc";
    std::string filename_2 = "temp_scenarioData_2.nc";

    simpleNetCDF(filename, true);
    simpleNetCDF(filename_2, false);

    try {
        NcFile dataFile(filename, NcFile::read);
        NcFile dataFile_2(filename_2, NcFile::read);
        DiscreteScenarioSet dss;
        DiscreteScenarioSet dss_2;

        dss.deserialize(dataFile);
        dss_2.deserialize(dataFile_2);

        std::cout << "    Simple checks" << std::endl;
        assert( dss.nbScenarios == 5 );
        assert( dss.scenarioSize == 10 ); 

        std::cout << "    Default probability distribution is uniform" << std::endl;
        assert( dss_2.scenarioProbabilities.size() == 5);
        assert( dss_2.scenarioProbabilities[3] = 0.2);

        // Delete the files after testing
        remove(filename.c_str());
        remove(filename_2.c_str());

        std::cout << "    Deserialization test passed" << std::endl;

    } catch (const std::exception& e) {
        // Ensure the file is deleted even if an exception occurs
        remove(filename.c_str());
        remove(filename_2.c_str());
        std::cerr << "    An error occurred: " << e.what() << std::endl;
    }
}

/*--------------------------------------------------------------------------*/
/*------------------------------- MAIN -------------------------------------*/
/*--------------------------------------------------------------------------*/

int main(int argc, char** argv) {

    // In the absence of the -verbose command, suppress all std::cout
    bool verbose = false; 
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "-verbose") {
            verbose = true;
        }
    }

    // Scope in which cout is conditionally allowed by the verbose variable
    { 
      CoutSuppressor out(verbose);

      // Unit tests
      test_deserialize();
    }

    // This message will always print
    std::cout << "Test_discret done" << std::endl;

    return 0;
}
