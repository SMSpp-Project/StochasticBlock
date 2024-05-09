#include <iostream>
#include <netcdf>
#include "DiscreteScenarioSet.h"

using namespace netCDF;
using namespace SMSpp_di_unipi_it;

// Function to create a netCDF file
void createNetCDF(const std::string& filename) {
    NcFile dataFile(filename, NcFile::replace);
    NcDim scenarioDim = dataFile.addDim("scenario", 5);
    NcDim dimensionDim = dataFile.addDim("dimension", 10);
    NcVar dataVar = dataFile.addVar("scenarios", ncDouble, {scenarioDim, dimensionDim});
    NcVar probVar = dataFile.addVar("probabilities", ncDouble, scenarioDim);

    std::vector<double> scenarios(scenarioDim.getSize() * dimensionDim.getSize(), 0.5);
    std::vector<double> probabilities(scenarioDim.getSize(), 1.0 / scenarioDim.getSize());

    dataVar.putVar(scenarios.data());
    probVar.putVar(probabilities.data());
}

int main() {
    std::string filename = "temp_scenarioData.nc";
    createNetCDF(filename);

    try {
        NcFile dataFile(filename, NcFile::read);
        DiscreteScenarioSet dss;
        dss.deserialize(dataFile);

        // Will add more tests
        std::cout << "Deserialization successful." << std::endl;

        // Delete the file after testing
        remove(filename.c_str());
        std::cout << "Temporary file deleted." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "An error occurred: " << e.what() << std::endl;
        // Ensure the file is deleted even if an exception occurs
        remove(filename.c_str());
    }

    return 0;
}
