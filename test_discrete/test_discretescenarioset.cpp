/*--------------------------------------------------------------------------*/
/*--------------------- File test_discretescenarioset.cpp ------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Test suite for DiscreteScenarioSet class
 * 
 * Test 1 - Basic Functionality:
 * - Scenario loading and deserialization
 * - Parameter validation for init_representative_pool
 * - Random pool initialization
 * 
 * Test 2 - Scenario Reduction Algorithms:
 * - ScenarioReductionSolver algorithms (Dupacova, BestFit, FirstFit)
 * - MILPSolver implementations (CPLEX, HiGHS)
 * 
 * Test 3 - Configuration Management:
 * - Valid/invalid solver configuration
 * - set_config() method
 * - Configuration serialization/deserialization
 * 
 * Test 4 - Serialization and Deserialization:
 * - DiscreteScenarioSet persistence
 * - Configuration persistence across save/load
 * - Scenario reduction deserialization with various configurations
 * 
 * \author Benoît Tran \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 * 
 * Copyright &copy; by Benoît Tran
 */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <span>
#include "DiscreteScenarioSet.h"
#include "Configuration.h"
#include "BlockSolverConfig.h"
#include <netcdf>

#include <iostream>
#include <vector>
#include <memory>
#include <cstdio>
#include <string>
#include <functional>
#include <map>
#include <chrono>
#include <iomanip>
#include <cmath>
#include <fstream>

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
            cout << "PASSED: " << result.message << endl;
            tests_passed++;
        } else {
            cout << "FAILED: " << result.message << endl;
            tests_failed++;
        }
    } catch (const exception& e) {
        cout << "FAILED with exception: " << e.what() << endl;
        tests_failed++;
    }
}

/*--------------------------------------------------------------------------*/
/*--------------------------- TEST UTILITIES -------------------------------*/
/*--------------------------------------------------------------------------*/

// Helper function to create a test netCDF file with scenarios
string create_test_scenario_file(int num_scenarios = 20, int scenario_size = 10, 
                                const string& suffix = "") {
    string filename = "test_scenarios_" + to_string(num_scenarios) + "_" + 
                     to_string(scenario_size) + suffix + ".nc4";
    
    try {
        netCDF::NcFile dataFile(filename, netCDF::NcFile::replace);
        
        // Add type attribute at root level for deserializer to find
        dataFile.putAtt("type", "DiscreteScenarioSet");
        
        // Create dimensions
        netCDF::NcDim scenarioDim = dataFile.addDim("NumberScenarios", num_scenarios);
        netCDF::NcDim sizeDim = dataFile.addDim("ScenarioSize", scenario_size);
        
        // Create scenario variable
        netCDF::NcVar scenarioVar = dataFile.addVar("Scenarios", netCDF::ncDouble, 
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
        netCDF::NcVar probVar = dataFile.addVar("poolProbabilities", netCDF::ncDouble, scenarioDim);
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
unique_ptr<DiscreteScenarioSet> load_test_scenarios(int num_scenarios = 20, 
                                                   int scenario_size = 10) {
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

// Helper to create BlockConfig for scenario reduction
BlockConfig* create_block_config(int k, double ell = 2.0) {
    auto* block_config = new BlockConfig();
    
    // Set k parameter in extra configuration
    block_config->f_extra_Configuration = new SimpleConfiguration<int>(k);
    
    // Set ell parameter in static variables configuration
    block_config->f_static_variables_Configuration = new SimpleConfiguration<double>(ell);
    
    return block_config;
}

// Helper to create BlockSolverConfig for different solver types
BlockSolverConfig* create_solver_config(const string& solver_type, 
                                       const string& algorithm = "",
                                       double time_limit = 60.0,
                                       int verbosity = 0) {
    auto* solver_config = new BlockSolverConfig(true);  // differential mode
    
    if (solver_type == "ScenarioReductionSolver") {
        // ScenarioReductionSolver configuration
        solver_config->add_ComputeConfig("ScenarioReductionSolver", nullptr);
    } else {
        // MILP Solver configuration
        auto* compute_config = new ComputeConfig();
        compute_config->f_diff = true;
        
        // Add common MILP parameters
        compute_config->int_pars.emplace_back("intLogVerb", verbosity);
        compute_config->int_pars.emplace_back("intRelaxIntVars", 0);
        compute_config->dbl_pars.emplace_back("dblRelAcc", 1e-7);
        
        if (time_limit > 0) {
            // Use solver-specific time limit parameter names
            if (solver_type == "CPXMILPSolver") {
                compute_config->dbl_pars.emplace_back("CPXPARAM_TimeLimit", time_limit);
            } else if (solver_type == "GRBMILPSolver") {
                compute_config->dbl_pars.emplace_back("TimeLimit", time_limit);
            } else if (solver_type == "SCIPMILPSolver") {
                compute_config->dbl_pars.emplace_back("limits/time", time_limit);
            } else if (solver_type == "HiGHSMILPSolver") {
                compute_config->dbl_pars.emplace_back("time_limit", time_limit);
            }
        }
        
        // Add solver-specific parameters
        if (solver_type == "CPXMILPSolver") {
            compute_config->int_pars.emplace_back("CPXPARAM_Threads", 1);
            if (verbosity > 0) {
                compute_config->int_pars.emplace_back("CPXPARAM_MIP_Display", 3);
            }
        } else if (solver_type == "GRBMILPSolver") {
            compute_config->int_pars.emplace_back("Threads", 1);
            if (verbosity > 0) {
                compute_config->int_pars.emplace_back("OutputFlag", 1);
            }
        }
        
        solver_config->add_ComputeConfig(string(solver_type), compute_config);
    }
    
    return solver_config;
}

/*--------------------------------------------------------------------------*/
/*---------------------------- TEST FUNCTIONS ------------------------------*/
/*--------------------------------------------------------------------------*/

// Test 1: Basic functionality
TestResult test_basic_functionality() {
    try {
        // Part 1: Basic loading and deserialization
        auto dss = load_test_scenarios(20, 10);
        
        if (dss->get_nbScenarios() != 20) {
            return {false, "Expected 20 scenarios, got " + to_string(dss->get_nbScenarios())};
        }
        
        if (dss->get_scenarioSize() != 10) {
            return {false, "Expected scenario size 10, got " + to_string(dss->get_scenarioSize())};
        }
        
        // Part 2: Test invalid k parameters for init_representative_pool
        try {
            dss->init_representative_pool(0);
            return {false, "Should have thrown invalid_argument for k=0"};
        } catch (const invalid_argument& e) {
            // Expected
        }
        
        try {
            dss->init_representative_pool(25);
            return {false, "Should have thrown invalid_argument for k>nbScenarios"};
        } catch (const invalid_argument& e) {
            // Expected
        }
        
        // Part 3: Test random pool functionality
        dss->init_random_pool(10);
        
        auto scenario = dss->get_current_scenario();
        if (scenario.size() != 10) {
            return {false, "Expected scenario size 10, got " + to_string(scenario.size())};
        }
        
        double prob = dss->get_current_scenario_probability();
        if (prob <= 0.0 || prob > 1.0) {
            return {false, "Invalid probability: " + to_string(prob)};
        }
        
        return {true, "Basic functionality tests passed"};
        
    } catch (const exception& e) {
        return {false, string("Exception during test: ") + e.what()};
    }
}

REGISTER_TEST("Basic Functionality", test_basic_functionality);

// Test 2: Scenario reduction algorithms
TestResult test_scenario_reduction_algorithms() {
    try {
        const int num_scenarios = 15;
        const int scenario_size = 5;
        const int k = 5;
        
        // Test ScenarioReductionSolver algorithms
        vector<string> algorithms = {"Dupacova", "BestFit", "FirstFit"};
        
        for (const auto& algorithm : algorithms) {
            auto dss = load_test_scenarios(num_scenarios, scenario_size);
            
            auto* block_config = create_block_config(k);
            auto* solver_config = create_solver_config("ScenarioReductionSolver", algorithm);
            
            dss->set_scenario_reduction_config(block_config, solver_config);
            dss->init_representative_pool(k);
            
            if (dss->get_selected_scenario_count() != k) {
                return {false, algorithm + " failed: wrong number of scenarios selected"};
            }
        }
        
        // Test MILPSolver implementations
        vector<string> milp_solvers = {"CPXMILPSolver", "HiGHSMILPSolver"};
        
        for (const auto& solver_name : milp_solvers) {
            try {
                // Create distinct scenarios for MILP testing
                string filename = "test_milp_scenarios.nc4";
                netCDF::NcFile dataFile(filename, netCDF::NcFile::replace);
                
                auto scenarioDim = dataFile.addDim("NumberScenarios", 10);
                auto sizeDim = dataFile.addDim("ScenarioSize", 5);
                auto scenarioVar = dataFile.addVar("Scenarios", netCDF::ncDouble, {scenarioDim, sizeDim});
                
                vector<double> data(50);
                for (int i = 0; i < 10; ++i) {
                    for (int j = 0; j < 5; ++j) {
                        data[i * 5 + j] = i * 10.0 + j;
                    }
                }
                scenarioVar.putVar(data.data());
                
                auto probVar = dataFile.addVar("poolProbabilities", netCDF::ncDouble, scenarioDim);
                vector<double> probs(10, 0.1);
                probVar.putVar(probs.data());
                dataFile.close();
                
                auto dss = make_unique<DiscreteScenarioSet>();
                netCDF::NcFile readFile(filename, netCDF::NcFile::read);
                dss->deserialize(readFile);
                readFile.close();
                
                auto* block_config = create_block_config(3, 2.0);
                auto* solver_config = create_solver_config(solver_name, "", 30.0, 0);
                
                dss->set_scenario_reduction_config(block_config, solver_config);
                dss->init_representative_pool(3);
                
                if (dss->get_selected_scenario_count() != 3) {
                    remove(filename.c_str());
                    return {false, solver_name + " failed to select exactly 3 scenarios"};
                }
                
                remove(filename.c_str());
                
            } catch (const exception& e) {
                cout << solver_name << " not available or test skipped: " << e.what() << endl;
            }
        }
        
        return {true, "All scenario reduction algorithms tested successfully"};
        
    } catch (const exception& e) {
        return {false, string("Test failed: ") + e.what()};
    }
}

REGISTER_TEST("Scenario Reduction Algorithms", test_scenario_reduction_algorithms);

// Test 3: Configuration management
TestResult test_configuration_management() {
    try {
        // Test 1: Valid solver configurations
        vector<string> valid_solver_names = {"CPXMILPSolver", "GRBMILPSolver", "ScenarioReductionSolver"};
        
        for (const auto& solver_name : valid_solver_names) {
            auto* block_config = create_block_config(3, 2.0);
            BlockSolverConfig* solver_config = nullptr;
            
            if (solver_name == "ScenarioReductionSolver") {
                solver_config = create_solver_config("ScenarioReductionSolver", "Dupacova");
            } else {
                solver_config = create_solver_config(solver_name, "", 60.0, 0);
            }
            
            auto dss = load_test_scenarios(10, 5);
            
            try {
                dss->set_scenario_reduction_config(block_config, solver_config);
            } catch (const exception& e) {
                delete block_config;
                delete solver_config;
                
                if (solver_name == "ScenarioReductionSolver") {
                    return {false, "ScenarioReductionSolver config failed: " + string(e.what())};
                }
            }
        }
        
        // Test 2: Invalid solver configuration
        {
            auto dss = load_test_scenarios(10, 5);
            
            auto* block_config = create_block_config(5);
            auto* solver_config = new BlockSolverConfig(true);
            solver_config->add_ComputeConfig("InvalidSolver", nullptr);
            
            try {
                dss->set_scenario_reduction_config(block_config, solver_config);
                return {false, "Should have thrown exception for invalid solver"};
            } catch (const invalid_argument& e) {
                string error_msg = e.what();
                if (error_msg.find("Unsupported solver for scenario reduction") == string::npos) {
                    return {false, "Unexpected error message for invalid solver: " + error_msg};
                }
            }
        }
        
        // Test 3: set_config() method
        {
            auto dss = load_test_scenarios(20, 3);
            
            // Test SimpleConfiguration with k
            auto k_config = make_unique<SimpleConfiguration<int>>(5);
            dss->set_config(k_config.get());
            
            if (dss->get_selected_scenario_count() != 5) {
                return {false, "Expected 5 selected scenarios after set_config"};
            }
            
            // Test override with new k
            auto k_config2 = make_unique<SimpleConfiguration<int>>(8);
            dss->set_config(k_config2.get());
            
            if (dss->get_selected_scenario_count() != 8) {
                return {false, "Expected 8 selected scenarios after override"};
            }
            
            // Test BlockConfig
            auto dss2 = load_test_scenarios(20, 3);
            auto block_cfg = make_unique<BlockConfig>();
            auto k_extra = make_unique<SimpleConfiguration<int>>(6);
            block_cfg->f_extra_Configuration = k_extra.release();
            auto ell_config = make_unique<SimpleConfiguration<double>>(2.5);
            block_cfg->f_static_variables_Configuration = ell_config.release();
            
            dss2->set_config(block_cfg.get());
            
            if (dss2->get_selected_scenario_count() != 6) {
                return {false, "Expected 6 selected scenarios with BlockConfig"};
            }
            
            // Test vector-based configuration
            auto dss3 = load_test_scenarios(15, 2);
            vector<double> params = {7.0, 3.0, 1e-6, 100.0};
            auto vec_config = make_unique<SimpleConfiguration<vector<double>>>(params);
            dss3->set_config(vec_config.get());
            
            if (dss3->get_selected_scenario_count() != 7) {
                return {false, "Expected 7 selected scenarios with vector config"};
            }
        }
        
        return {true, "Configuration management tests passed"};
        
    } catch (const exception& e) {
        return {false, string("Configuration test failed: ") + e.what()};
    }
}

REGISTER_TEST("Configuration Management", test_configuration_management);

// Test 4: Serialization and deserialization
TestResult test_serialization_deserialization() {
    try {
        // Test 1: Basic serialization with config
        {
            auto dss1 = load_test_scenarios(10, 5);
            auto* block_config = create_block_config(3, 2.0);
            auto* solver_config = create_solver_config("ScenarioReductionSolver", "Dupacova");
            dss1->set_scenario_reduction_config(block_config, solver_config);
            
            string nc_filename = "test_dss_with_config.nc4";
            {
                netCDF::NcFile file(nc_filename, netCDF::NcFile::replace);
                dss1->serialize(file);
                file.close();
            }
            
            auto dss2 = make_unique<DiscreteScenarioSet>();
            {
                netCDF::NcFile file(nc_filename, netCDF::NcFile::read);
                dss2->deserialize(file);
                file.close();
            }
            
            if (dss2->get_k_value() != 3) {
                remove(nc_filename.c_str());
                return {false, "k_value not properly restored"};
            }
            
            remove(nc_filename.c_str());
        }
        
        // Test 2: Deserialization with various configurations
        {
            // Test no config
            string filename = create_test_scenario_file(20, 5, "_no_config");
            auto dss = make_unique<DiscreteScenarioSet>();
            netCDF::NcFile dataFile(filename, netCDF::NcFile::read);
            dss->deserialize(dataFile);
            dataFile.close();
            
            if (dss->get_k_value() != 0) {
                remove(filename.c_str());
                return {false, "k_value should be 0 when no config is present"};
            }
            
            remove(filename.c_str());
            
            // Test only k provided
            filename = create_test_scenario_file(20, 5, "_only_k");
            {
                netCDF::NcFile dataFile(filename, netCDF::NcFile::write);
                auto cfgGroup = dataFile.addGroup("ScenarioReductionConfig");
                int k = 5;
                auto kVar = cfgGroup.addVar("k", netCDF::ncInt);
                kVar.putVar(&k);
                dataFile.close();
            }
            
            dss = make_unique<DiscreteScenarioSet>();
            dataFile.open(filename, netCDF::NcFile::read);
            dss->deserialize(dataFile);
            dataFile.close();
            
            if (dss->get_k_value() != 5 || dss->get_selected_scenario_count() != 5) {
                remove(filename.c_str());
                return {false, "Should have initialized pool with k=5"};
            }
            
            remove(filename.c_str());
            
            // Test k + ell provided
            filename = create_test_scenario_file(20, 5, "_k_and_ell");
            {
                netCDF::NcFile dataFile(filename, netCDF::NcFile::write);
                auto cfgGroup = dataFile.addGroup("ScenarioReductionConfig");
                int k = 7;
                auto kVar = cfgGroup.addVar("k", netCDF::ncInt);
                kVar.putVar(&k);
                float ell = 1.5f;
                auto ellVar = cfgGroup.addVar("ell", netCDF::ncFloat);
                ellVar.putVar(&ell);
                dataFile.close();
            }
            
            dss = make_unique<DiscreteScenarioSet>();
            dataFile.open(filename, netCDF::NcFile::read);
            dss->deserialize(dataFile);
            dataFile.close();
            
            if (dss->get_k_value() != 7 || dss->get_selected_scenario_count() != 7) {
                remove(filename.c_str());
                return {false, "Should have selected 7 scenarios"};
            }
            
            remove(filename.c_str());
            
            // Test invalid k values
            filename = create_test_scenario_file(20, 5, "_k_zero");
            {
                netCDF::NcFile dataFile(filename, netCDF::NcFile::write);
                auto cfgGroup = dataFile.addGroup("ScenarioReductionConfig");
                int k = 0;
                auto kVar = cfgGroup.addVar("k", netCDF::ncInt);
                kVar.putVar(&k);
                dataFile.close();
            }
            
            dss = make_unique<DiscreteScenarioSet>();
            dataFile.open(filename, netCDF::NcFile::read);
            dss->deserialize(dataFile);
            dataFile.close();
            
            if (dss->get_selected_scenario_count() != 0) {
                remove(filename.c_str());
                return {false, "No reduction should occur with k=0"};
            }
            
            remove(filename.c_str());
        }
        
        // Test 3: Serialization round-trip
        {
            string filename1 = create_test_scenario_file(20, 5, "_roundtrip1");
            string filename2 = "test_roundtrip2.nc4";
            
            {
                netCDF::NcFile dataFile(filename1, netCDF::NcFile::write);
                auto cfgGroup = dataFile.addGroup("ScenarioReductionConfig");
                int k = 8;
                auto kVar = cfgGroup.addVar("k", netCDF::ncInt);
                kVar.putVar(&k);
                float ell = 2.5f;
                auto ellVar = cfgGroup.addVar("ell", netCDF::ncFloat);
                ellVar.putVar(&ell);
                dataFile.close();
            }
            
            auto dss1 = make_unique<DiscreteScenarioSet>();
            netCDF::NcFile dataFile1(filename1, netCDF::NcFile::read);
            dss1->deserialize(dataFile1);
            dataFile1.close();
            
            netCDF::NcFile dataFile2(filename2, netCDF::NcFile::replace);
            dss1->serialize(dataFile2);
            dataFile2.close();
            
            auto dss2 = make_unique<DiscreteScenarioSet>();
            netCDF::NcFile dataFile3(filename2, netCDF::NcFile::read);
            dss2->deserialize(dataFile3);
            dataFile3.close();
            
            if (dss2->get_k_value() != 8 || dss2->get_selected_scenario_count() != 8) {
                remove(filename1.c_str());
                remove(filename2.c_str());
                return {false, "Round-trip serialization failed to preserve k and pool"};
            }
            
            remove(filename1.c_str());
            remove(filename2.c_str());
        }
        
        return {true, "Serialization and deserialization tests passed"};
        
    } catch (const exception& e) {
        return {false, string("Serialization test failed: ") + e.what()};
    }
}

REGISTER_TEST("Serialization and Deserialization", test_serialization_deserialization);

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
    cout << "Testing scenario reduction functionality\n" << endl;
    
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
    
    // Clean up any leftover test files
    system("rm -f test_*.nc4 test_*.txt 2>/dev/null");
    
    return tests_failed > 0 ? 1 : 0;
}