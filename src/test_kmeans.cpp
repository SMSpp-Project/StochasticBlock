#include <boost/multi_array.hpp>
#include <Eigen/Dense>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include "Kmeans.cpp" 

void testKmeans() {
    // Setting up the test data
    int numScenarios = 10;  // Number of data points
    int dimension = 3;      // Dimension of each data point
    boost::multi_array<double, 2> scenarios(boost::extents[numScenarios][dimension]);

    // Filling the test data with random values
    std::srand(static_cast<unsigned int>(std::time(nullptr)));
    for (int i = 0; i < numScenarios; ++i) {
        for (int j = 0; j < dimension; ++j) {
            scenarios[i][j] = -2.5+5*static_cast<double>(rand()) / RAND_MAX;  // Random double between -2.5 and 2.5
        }
    }

    // Instantiating Kmeans
    SMSpp_di_unipi_it::Kmeans kmeans_instance(scenarios);

    // Run kMeans clustering
    int k = 3;  // Number of clusters
    auto [labels, centers] = kmeans_instance.kMeans(k, kmeans_instance.poolMap);  

    // Output the results
    std::cout << "Cluster labels:\n";
    for (int i = 0; i < numScenarios; ++i) {
        std::cout << "Scenario " << i << ": Cluster " << labels[i] << std::endl;
    }

    // Output the centers
    for (const auto& center : centers) {
        std::cout << "Center: " << center.transpose() << std::endl;
    }
}

int main() {
    testKmeans();
    return 0;
}


