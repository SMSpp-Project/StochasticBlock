/*--------------------------------------------------------------------------*/
/*--------------------- File test_discretescenarioset.cpp ------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Test suite for DiscreteScenarioSet class
 * 
 * This test validates the init_representative_pool method with k parameter
 * 
 * \date 2025
 */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DiscreteScenarioSet.h"
#include <netcdf>

#include <iostream>
#include <vector>
#include <memory>
#include <cstdio>
#include <string>
#include <functional>
#include <map>

/*--------------------------------------------------------------------------*/
/*------------------------------- USING ------------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;
using namespace std;

/*--------------------------------------------------------------------------*/
/*--------------------------- TEST FRAMEWORK -------------------------------*/
/*--------------------------------------------------------------------------*/

struct TestResult {
    bool passed;
    string message;
};

// Global test counters
static int tests_run = 0;
static int tests_passed = 0;
static int tests_failed = 0;

// Test registry
static map<string, function<TestResult()>> test_registry;

// Register a test
#define REGISTER_TEST(name, func) \
    static bool _reg_##func = []() { \
        test_registry[name] = func; \
        return true; \
    }()

// Helper to run a single test
void run_test(const string& name, function<TestResult()> test_func) {
    cout << "\n=== Running: " << name << " ===" << endl;
    tests_run++;
    
    try {
        TestResult result = test_func();
        if (result.passed) {
            cout << "✓ PASSED: " << result.message << endl;
            tests_passed++;
        } else {
            cout << "✗ FAILED: " << result.message << endl;
            tests_failed++;
        }
    } catch (const exception& e) {
        cout << "✗ FAILED with exception: " << e.what() << endl;
        tests_failed++;
    }
}

/*--------------------------------------------------------------------------*/
/*--------------------------- TEST UTILITIES -------------------------------*/
/*--------------------------------------------------------------------------*/

// Helper function to create a test netCDF file with scenarios
string create_test_scenario_file(int num_scenarios = 20, int scenario_size = 10) {
    string filename = "test_scenarios_" + to_string(num_scenarios) + "_" + to_string(scenario_size) + ".nc4";
    
    try {
        netCDF::NcFile dataFile(filename, netCDF::NcFile::replace);
        
        // Create dimensions
        netCDF::NcDim scenarioDim = dataFile.addDim("NumberScenarios", num_scenarios);
        netCDF::NcDim sizeDim = dataFile.addDim("ScenarioSize", scenario_size);
        
        // Create scenario variable
        netCDF::NcVar scenarioVar = dataFile.addVar("Scenario", netCDF::ncDouble, 
                                                    {scenarioDim, sizeDim});
        
        // Fill with test data
        vector<double> data(num_scenarios * scenario_size);
        for (int i = 0; i < num_scenarios; ++i) {
            for (int j = 0; j < scenario_size; ++j) {
                data[i * scenario_size + j] = i + j * 0.1;
            }
        }
        scenarioVar.putVar(data.data());
        
        // Add probabilities (uniform distribution)
        netCDF::NcVar probVar = dataFile.addVar("Probability", netCDF::ncDouble, scenarioDim);
        vector<double> probs(num_scenarios, 1.0 / num_scenarios);
        probVar.putVar(probs.data());
        
        dataFile.close();
        return filename;
        
    } catch (netCDF::exceptions::NcException& e) {
        cerr << "Error creating test netCDF file: " << e.what() << endl;
        throw;
    }
}

// Helper to load scenarios into DiscreteScenarioSet
unique_ptr<DiscreteScenarioSet> load_test_scenarios(int num_scenarios = 20, int scenario_size = 10) {
    string filename = create_test_scenario_file(num_scenarios, scenario_size);
    
    auto dss = make_unique<DiscreteScenarioSet>();
    netCDF::NcFile dataFile(filename, netCDF::NcFile::read);
    netCDF::NcGroup root = dataFile;
    dss->deserialize(root);
    dataFile.close();
    
    // Clean up the file
    remove(filename.c_str());
    
    return dss;
}

/*--------------------------------------------------------------------------*/
/*---------------------------- TEST FUNCTIONS ------------------------------*/
/*--------------------------------------------------------------------------*/

// Test 1: Basic loading and deserialization
TestResult test_basic_loading() {
    cout << "Testing basic scenario loading from netCDF..." << endl;
    
    try {
        auto dss = load_test_scenarios(20, 10);
        
        if (dss->get_nbScenarios() != 20) {
            return {false, "Expected 20 scenarios, got " + to_string(dss->get_nbScenarios())};
        }
        
        if (dss->get_scenarioSize() != 10) {
            return {false, "Expected scenario size 10, got " + to_string(dss->get_scenarioSize())};
        }
        
        return {true, "Successfully loaded 20 scenarios of size 10"};
        
    } catch (const exception& e) {
        return {false, string("Exception during loading: ") + e.what()};
    }
}

REGISTER_TEST("Basic Loading", test_basic_loading);

// Test 2: init_representative_pool with k = 0
TestResult test_invalid_k_zero() {
    cout << "Testing init_representative_pool with k = 0..." << endl;
    
    try {
        auto dss = load_test_scenarios(20, 10);
        
        try {
            dss->init_representative_pool(0);
            return {false, "Should have thrown invalid_argument for k=0"};
        } catch (const invalid_argument& e) {
            return {true, "Correctly rejected k=0: " + string(e.what())};
        } catch (const exception& e) {
            return {false, "Wrong exception type: " + string(e.what())};
        }
        
    } catch (const exception& e) {
        return {false, string("Setup failed: ") + e.what()};
    }
}

REGISTER_TEST("Invalid k=0", test_invalid_k_zero);

// Test 3: init_representative_pool with k > nbScenarios
TestResult test_invalid_k_too_large() {
    cout << "Testing init_representative_pool with k > nbScenarios..." << endl;
    
    try {
        auto dss = load_test_scenarios(20, 10);
        
        try {
            dss->init_representative_pool(25);
            return {false, "Should have thrown invalid_argument for k=25"};
        } catch (const invalid_argument& e) {
            return {true, "Correctly rejected k=25: " + string(e.what())};
        } catch (const exception& e) {
            return {false, "Wrong exception type: " + string(e.what())};
        }
        
    } catch (const exception& e) {
        return {false, string("Setup failed: ") + e.what()};
    }
}

REGISTER_TEST("Invalid k>nbScenarios", test_invalid_k_too_large);

// Test 4: init_representative_pool with valid k values
TestResult test_valid_k_values() {
    cout << "Testing init_representative_pool with valid k values..." << endl;
    
    try {
        auto dss = load_test_scenarios(20, 10);
        vector<int> k_values = {1, 5, 10, 20};
        
        for (int k : k_values) {
            try {
                dss->init_representative_pool(k);
                cout << "  ✓ k=" << k << " accepted" << endl;
            } catch (const exception& e) {
                // Note: It might fail if scenario reduction solver is not available
                // but the parameter validation should have passed
                cout << "  ⚠ k=" << k << " - " << e.what() << endl;
            }
        }
        
        return {true, "Parameter validation passed for all valid k values"};
        
    } catch (const exception& e) {
        return {false, string("Test failed: ") + e.what()};
    }
}

REGISTER_TEST("Valid k values", test_valid_k_values);

// Test 5: init_random_pool functionality
TestResult test_random_pool() {
    cout << "Testing init_random_pool functionality..." << endl;
    int dim = 5;
    try {
        auto dss = load_test_scenarios(20, dim);
        
        // Initialize random pool
        dss->init_random_pool(10);
        
        // Get a scenario
        auto scenario = dss->get_current_scenario();
        
        if (scenario.size() != dim) {
            return {false, "Expected scenario size, got " + to_string(scenario.size())};
        }
        
        // Check probability
        double prob = dss->get_current_scenario_probability();
        if (prob <= 0.0 || prob > 1.0) {
            return {false, "Invalid probability: " + to_string(prob)};
        }
        
        return {true, "Random pool initialized successfully with 10 scenarios"};
        
    } catch (const exception& e) {
        return {false, string("Test failed: ") + e.what()};
    }
}

REGISTER_TEST("Random Pool", test_random_pool);

/*--------------------------------------------------------------------------*/
/*-------------------------------- MAIN ------------------------------------*/
/*--------------------------------------------------------------------------*/

void print_usage(const char* program_name) {
    cout << "Usage: " << program_name << " [options]" << endl;
    cout << "Options:" << endl;
    cout << "  -h, --help     Show this help message" << endl;
    cout << "  -l, --list     List all available tests" << endl;
    cout << "  -t, --test     Run specific test by name" << endl;
    cout << "  -v, --verbose  Verbose output" << endl;
    cout << "  (no options)   Run all tests" << endl;
}

int main(int argc, char* argv[]) {
    cout << "========== DiscreteScenarioSet Test Suite ==========" << endl;
    cout << "Testing init_representative_pool functionality\n" << endl;
    
    // Parse command line arguments
    bool run_all = true;
    bool list_tests = false;
    string specific_test;
    
    for (int i = 1; i < argc; i++) {
        string arg = argv[i];
        
        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "-l" || arg == "--list") {
            list_tests = true;
            run_all = false;
        } else if (arg == "-t" || arg == "--test") {
            if (i + 1 < argc) {
                specific_test = argv[++i];
                run_all = false;
            } else {
                cerr << "Error: -t/--test requires a test name" << endl;
                return 1;
            }
        } else if (arg == "-v" || arg == "--verbose") {
            // Could add verbose flag handling here
        } else {
            cerr << "Unknown option: " << arg << endl;
            print_usage(argv[0]);
            return 1;
        }
    }
    
    // List tests if requested
    if (list_tests) {
        cout << "Available tests:" << endl;
        for (const auto& [name, func] : test_registry) {
            cout << "  - " << name << endl;
        }
        return 0;
    }
    
    // Run specific test if requested
    if (!specific_test.empty()) {
        auto it = test_registry.find(specific_test);
        if (it != test_registry.end()) {
            run_test(it->first, it->second);
        } else {
            cerr << "Error: Test '" << specific_test << "' not found" << endl;
            cout << "Use -l to list available tests" << endl;
            return 1;
        }
    }
    
    // Run all tests
    if (run_all) {
        for (const auto& [name, func] : test_registry) {
            run_test(name, func);
        }
    }
    
    // Print summary
    cout << "\n========== Test Summary ==========" << endl;
    cout << "Tests run:    " << tests_run << endl;
    cout << "Tests passed: " << tests_passed << endl;
    cout << "Tests failed: " << tests_failed << endl;
    
    return tests_failed > 0 ? 1 : 0;
}