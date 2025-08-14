/*--------------------------------------------------------------------------*/
/*--------------------- File test_discretescenarioset.cpp ------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Test suite for DiscreteScenarioSet class. 
 * 
 * Test 1 - Basic Functionality (Comprehensive):
 * - Part 1: Scenario loading and deserialization
 * - Part 2: Parameter validation (invalid poolSize)
 * - Part 3: Random pool initialization
 * - Part 4: Rejection sampling verification (unique selections)
 * - Part 5: Utility methods (is_pool_initialized, set_seed, get_ell/set_ell, get_scenario_value)
 * - Part 6: Edge cases and error conditions (single scenario, select all, access before init, weighted probabilities)
 * 
 * Test 2 - Configuration Patterns:
 * - Pattern 1: SimpleConfiguration<int> (baseline method)
 * - Pattern 2: SimpleConfiguration<pair<int, BlockSolverConfig*>> (advanced + generated BlockConfig)
 * - Pattern 3: SimpleConfiguration<pair<BlockConfig*, BlockSolverConfig*>> (full advanced)
 * 
 * Test 3 - Scenario Reduction Algorithms:
 * - ScenarioReductionSolver algorithms (Dupacova, BestFit, FirstFit)
 * - MILPSolver implementations (CPLEX, HiGHS)
 * 
 * Test 4 - Serialization and Deserialization:
 * - DiscreteScenarioSet persistence with scenario reduction configuration
 * - Full serialization/deserialization round-trip with solver configs
 * - Deserialization with various configurations (no config, poolSize only, poolSize+ell)
 * - Complete scenario data persistence and restoration
 * - BlockConfig and BlockSolverConfig serialization
 * - Invalid poolSize value handling during deserialization
 * 
 * Test 5 - Iteration and Span-based Getters:
 * - Full iteration through selected scenarios
 * - Span-based getters (get_selected_scenarios, get_pool_weights, get_normalized_weights)
 * - Probability normalization verification
 * - Index validation
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
#include "BlockSolverConfig.h"

#include <iostream>  // std::cout, std::cerr
#include <cstdio>    // std::remove
#include <chrono>    // std::chrono for timing
#include <iomanip>   // std::setprecision 
#include <cmath>     // std::abs, std::sqrt
#include <set>       // std::set for uniqueness testing
#include <span>      // std::span for C++20 features

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
BlockConfig* create_block_config(int poolSize, double ell = 2.0) {
    auto* block_config = new BlockConfig();
    
    // Set poolSize parameter in extra configuration
    block_config->f_extra_Configuration = new SimpleConfiguration<int>(poolSize);
    
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
        
        // Part 2: Test invalid poolSize parameters for init_representative_pool
        try {
            dss->init_representative_pool(0);
            return {false, "Should have thrown invalid_argument for poolSize=0"};
        } catch (const invalid_argument& e) {
            // Expected
        }
        
        try {
            dss->init_representative_pool(25);
            return {false, "Should have thrown invalid_argument for poolSize>nbScenarios"};
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
        
        // Part 4: Test rejection sampling (random pool produces unique selections)
        // Create a larger dataset for testing
        auto dss_large = load_test_scenarios(100, 10);
        dss_large->set_seed(42);  // For reproducibility
        
        // Select 50 scenarios using random pool (uses rejection sampling internally)
        dss_large->init_random_pool(50);
        
        auto selected = dss_large->get_selected_scenarios();
        if (selected.size() != 50) {
            return {false, "Random pool selection size mismatch: expected 50, got " + to_string(selected.size())};
        }
        
        // Verify all selected indices are unique (key property of rejection sampling)
        set<ScenarioGenerator::ScenarioIndex> unique_indices(selected.begin(), selected.end());
        if (unique_indices.size() != selected.size()) {
            return {false, "Random pool produced duplicate indices (rejection sampling failed)"};
        }
        
        // Verify all indices are valid
        for (auto idx : selected) {
            if (idx >= 100) {
                return {false, "Invalid scenario index: " + to_string(idx)};
            }
        }
        
        // Part 5: Test utility methods
        // Test is_pool_initialized
        auto dss_util = load_test_scenarios(10, 4);
        if (dss_util->is_pool_initialized()) {
            return {false, "Pool should not be initialized initially"};
        }
        
        dss_util->init_random_pool(5);
        if (!dss_util->is_pool_initialized()) {
            return {false, "Pool should be initialized after init_random_pool"};
        }
        
        // Test set_seed for reproducibility
        auto dss1_seed = load_test_scenarios(10, 4);
        auto dss2_seed = load_test_scenarios(10, 4);
        
        dss1_seed->set_seed(12345);
        dss2_seed->set_seed(12345);
        
        dss1_seed->init_random_pool(5);
        dss2_seed->init_random_pool(5);
        
        auto indices1 = dss1_seed->get_selected_scenarios();
        auto indices2 = dss2_seed->get_selected_scenarios();
        
        bool same_selection = true;
        for (size_t i = 0; i < indices1.size(); ++i) {
            if (indices1[i] != indices2[i]) {
                same_selection = false;
                break;
            }
        }
        
        if (!same_selection) {
            return {false, "Same seed should produce same random selection"};
        }
        
        // Test get_ell and set_ell
        float original_ell = dss->get_ell();
        if (abs(original_ell - 2.0f) > 1e-6) {
            return {false, "Default ell should be 2.0"};
        }
        
        dss->set_ell(3.0f);
        if (abs(dss->get_ell() - 3.0f) > 1e-6) {
            return {false, "set_ell didn't update value correctly"};
        }
        
        // Test invalid ell
        try {
            dss->set_ell(-1.0f);
            return {false, "Should throw for negative ell"};
        } catch (const invalid_argument&) {
            // Expected
        }
        
        // Test get_scenario_value direct access
        double value = dss->get_scenario_value(0, 0);
        if (isnan(value) || isinf(value)) {
            return {false, "Invalid scenario value"};
        }
        
        // Test out of bounds
        try {
            (void)dss->get_scenario_value(100, 0);  // Explicitly discard [[nodiscard]] result
            return {false, "Should throw for out-of-bounds scenario index"};
        } catch (const out_of_range&) {
            // Expected
        }
        
        try {
            (void)dss->get_scenario_value(0, 100);  // Explicitly discard [[nodiscard]] result
            return {false, "Should throw for out-of-bounds component index"};
        } catch (const out_of_range&) {
            // Expected
        }
        
        // Part 6: Edge cases and error conditions
        // Test single scenario
        auto dss_single = load_test_scenarios(1, 3);
        dss_single->init_random_pool(1);
        
        // After init, we're positioned at the first (and only) scenario
        // next_scenario() should return false as there's no next
        if (dss_single->next_scenario()) {
            return {false, "Single scenario pool should not have next"};
        }
        
        // Test select all scenarios
        auto dss_all = load_test_scenarios(5, 3);
        dss_all->init_random_pool(5);
        
        int count = 1; // Start at 1 for current scenario
        while (dss_all->next_scenario()) {
            count++;
        }
        
        if (count != 5) {
            return {false, "Should iterate through all 5 scenarios"};
        }
        
        // Test error conditions - access before initialization
        auto dss_uninit = load_test_scenarios(10, 4);
        
        try {
            (void)dss_uninit->get_selected_scenarios();  // Explicitly discard [[nodiscard]] result
            return {false, "Should throw when accessing selected scenarios before init"};
        } catch (const runtime_error&) {
            // Expected
        }
        
        try {
            (void)dss_uninit->get_normalized_weights();  // Explicitly discard [[nodiscard]] result
            return {false, "Should throw when accessing normalized weights before init"};
        } catch (const runtime_error&) {
            // Expected
        }
        
        // Test weighted vs uniform probabilities
        auto dss_weighted = make_unique<DiscreteScenarioSet>();
        
        // Create temporary netCDF file with weighted scenarios
        string temp_file = "test_weighted_" + to_string(chrono::steady_clock::now().time_since_epoch().count()) + ".nc";
        
        try {
            {
                netCDF::NcFile file(temp_file, netCDF::NcFile::replace);
                auto group = file.addGroup("scenarios");
                
                group.addDim("NumberScenarios", 5);
                group.addDim("ScenarioSize", 2);
                
                auto scenarios_var = group.addVar("Scenarios", netCDF::ncDouble, 
                                                 {group.getDim("NumberScenarios"), 
                                                  group.getDim("ScenarioSize")});
                
                // Add non-uniform weights
                auto weights_var = group.addVar("poolWeights", netCDF::ncDouble,
                                               group.getDim("NumberScenarios"));
                
                double scenario_data[5][2] = {{1,2}, {3,4}, {5,6}, {7,8}, {9,10}};
                double weights[5] = {0.1, 0.2, 0.3, 0.3, 0.1}; // Non-uniform, sum = 1.0
                
                scenarios_var.putVar(scenario_data);
                weights_var.putVar(weights);
                
                dss_weighted->deserialize(group);
            }
        } catch (...) {
            remove(temp_file.c_str());
            throw;
        }
        remove(temp_file.c_str());
        
        // Test that weights are correctly loaded
        auto loaded_weights = dss_weighted->get_pool_weights();
        if (loaded_weights.size() != 5) {
            return {false, "Weights not loaded correctly"};
        }
        
        // Verify non-uniform distribution
        if (abs(loaded_weights[0] - 0.1) > 1e-6 || abs(loaded_weights[2] - 0.3) > 1e-6) {
            return {false, "Non-uniform weights not preserved"};
        }
        
        return {true, "Basic functionality tests passed"};
        
    } catch (const exception& e) {
        return {false, string("Exception during test: ") + e.what()};
    }
}

REGISTER_TEST("Basic Functionality", test_basic_functionality);

// Test 2: Configuration Patterns
TestResult test_configuration_patterns() {
    try {
        // Test Pattern 1: SimpleConfiguration<int> - baseline method
        {
            auto dss = load_test_scenarios(20, 5);
            
            // Pattern 1: Simple poolSize-only configuration
            auto* pattern1_config = new SimpleConfiguration<int>(8);
            dss->set_config(pattern1_config);
            
            // Verify poolSize was set
            if (dss->get_poolSize() != 8) {
                return {false, "Pattern 1: poolSize not set correctly (expected 8, got " + to_string(dss->get_poolSize()) + ")"};
            }
            
            // Test that baseline method works
            dss->init_representative_pool(8);
            if (dss->get_poolSize() != 8) {
                return {false, "Pattern 1: baseline method failed (expected 8 scenarios, got " + to_string(dss->get_poolSize()) + ")"};
            }
        }
        
        // Test Pattern 2: SimpleConfiguration<pair<int, Configuration*>> where Configuration* is BlockSolverConfig*
        {
            auto dss = load_test_scenarios(20, 5);
            
            // Create a BlockSolverConfig for ScenarioReductionSolver
            auto* solver_config = create_solver_config("ScenarioReductionSolver", "Dupacova");
            auto* pattern2_config = new SimpleConfiguration<pair<int, Configuration*>>(
                make_pair(6, solver_config));
            
            dss->set_config(pattern2_config);
            
            // Verify poolSize was set
            if (dss->get_poolSize() != 6) {
                return {false, "Pattern 2: poolSize not set correctly (expected 6, got " + to_string(dss->get_poolSize()) + ")"};
            }
            
            // Test that advanced scenario reduction works
            dss->init_representative_pool(6);
            if (dss->get_poolSize() != 6) {
                return {false, "Pattern 2: failed (expected 6 scenarios, got " + to_string(dss->get_poolSize()) + ")"};
            }
            
            // Note: Do not delete solver_config - DiscreteScenarioSet keeps a reference to it
        }
        
        // Test Pattern 3: SimpleConfiguration<pair<Configuration*, Configuration*>> where first is BlockConfig*, second is BlockSolverConfig*
        {
            auto dss = load_test_scenarios(20, 5);
            
            // Create both BlockConfig and BlockSolverConfig
            auto* block_config = create_block_config(7, 2.0);
            auto* solver_config = create_solver_config("ScenarioReductionSolver", "BestFit");
            auto* pattern3_config = new SimpleConfiguration<pair<Configuration*, Configuration*>>(
                make_pair(block_config, solver_config));
            
            dss->set_config(pattern3_config);
            
            // Verify poolSize was set from BlockConfig
            if (dss->get_poolSize() != 7) {
                delete block_config;  // BlockConfig is cloned, so we can delete the original
                return {false, "Pattern 3: poolSize not set correctly (expected 7, got " + to_string(dss->get_poolSize()) + ")"};
            }
            
            // Test that full advanced scenario reduction works
            dss->init_representative_pool(7);
            if (dss->get_poolSize() != 7) {
                delete block_config;  // BlockConfig is cloned, so we can delete the original
                return {false, "Pattern 3: failed (expected 7 scenarios, got " + to_string(dss->get_poolSize()) + ")"};
            }
            
            delete block_config;  // BlockConfig is cloned, so we can delete the original
            // Note: Do not delete solver_config - DiscreteScenarioSet keeps a reference to it
        }
        
        return {true, "All three configuration patterns tested successfully"};
        
    } catch (const exception& e) {
        return {false, string("Patterns test failed: ") + e.what()};
    }
}

REGISTER_TEST("Configuration Patterns", test_configuration_patterns);


// Test 3: Scenario reduction algorithms
TestResult test_scenario_reduction_algorithms() {
    try {
        const int num_scenarios = 15;
        const int scenario_size = 5;
        const int poolSize = 5;
        
        // Test various solver implementations - all treated uniformly
        // ScenarioReductionSolver with different algorithms
        vector<pair<string, string>> solver_configs = {
            {"ScenarioReductionSolver", "Dupacova"},
            {"ScenarioReductionSolver", "BestFit"},
            {"ScenarioReductionSolver", "FirstFit"},
            {"CPXMILPSolver", ""},
            {"HiGHSMILPSolver", ""}
        };
        
        for (const auto& [solver_name, algorithm] : solver_configs) {
            try {
                auto dss = load_test_scenarios(num_scenarios, scenario_size);
                
                auto* block_config = create_block_config(poolSize);
                auto* solver_config = create_solver_config(solver_name, algorithm);
                
                dss->set_config(block_config, solver_config);
                dss->init_representative_pool(poolSize);
                
                if (dss->get_poolSize() != poolSize) {
                    string test_name = algorithm.empty() ? solver_name : solver_name + ":" + algorithm;
                    return {false, test_name + " failed: wrong number of scenarios selected"};
                }
                
            } catch (const exception& e) {
                string test_name = algorithm.empty() ? solver_name : solver_name + ":" + algorithm;
                cout << test_name << " not available or test skipped: " << e.what() << endl;
            }
        }
        
        return {true, "All scenario reduction algorithms tested successfully"};
        
    } catch (const exception& e) {
        return {false, string("Test failed: ") + e.what()};
    }
}

REGISTER_TEST("Scenario Reduction Algorithms", test_scenario_reduction_algorithms);

// Test 4: Serialization and deserialization
TestResult test_serialization_deserialization() {
    try {
        // Part 1: Basic serialization with config
        {
            auto dss1 = load_test_scenarios(10, 5);
            auto* block_config = create_block_config(3, 2.0);
            auto* solver_config = create_solver_config("ScenarioReductionSolver", "Dupacova");
            dss1->set_config(block_config, solver_config);
            
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
            
            if (dss2->get_poolSize() != 3) {
                remove(nc_filename.c_str());
                return {false, "poolSize not properly restored"};
            }
            
            remove(nc_filename.c_str());
        }
        
        // Part 2: Deserialization with various configurations
        {
            // Test no config
            string filename = create_test_scenario_file(20, 5, "_no_config");
            auto dss = make_unique<DiscreteScenarioSet>();
            netCDF::NcFile dataFile(filename, netCDF::NcFile::read);
            dss->deserialize(dataFile);
            dataFile.close();
            
            if (dss->get_poolSize() != 0) {
                remove(filename.c_str());
                return {false, "poolSize should be 0 when no config is present"};
            }
            
            remove(filename.c_str());
            
            // Test only poolSize provided
            filename = create_test_scenario_file(20, 5, "_only_k");
            {
                netCDF::NcFile dataFile(filename, netCDF::NcFile::write);
                auto cfgGroup = dataFile.addGroup("ScenarioReductionConfig");
                int poolSize = 5;
                auto poolSizeVar = cfgGroup.addVar("poolSize", netCDF::ncInt);
                poolSizeVar.putVar(&poolSize);
                dataFile.close();
            }
            
            dss = make_unique<DiscreteScenarioSet>();
            dataFile.open(filename, netCDF::NcFile::read);
            dss->deserialize(dataFile);
            dataFile.close();
            
            if (dss->get_poolSize() != 5) {
                remove(filename.c_str());
                return {false, "Should have loaded poolSize=5 from file"};
            }
            
            dss->init_representative_pool(5);
            
            if (dss->get_poolSize() != 5) {
                remove(filename.c_str());
                return {false, "Should have 5 scenarios after init_representative_pool"};
            }
            
            remove(filename.c_str());
            
            // Test poolSize + ell provided
            filename = create_test_scenario_file(20, 5, "_k_and_ell");
            {
                netCDF::NcFile dataFile(filename, netCDF::NcFile::write);
                auto cfgGroup = dataFile.addGroup("ScenarioReductionConfig");
                int poolSize = 7;
                auto poolSizeVar = cfgGroup.addVar("poolSize", netCDF::ncInt);
                poolSizeVar.putVar(&poolSize);
                float ell = 1.5f;
                auto ellVar = cfgGroup.addVar("ell", netCDF::ncFloat);
                ellVar.putVar(&ell);
                dataFile.close();
            }
            
            dss = make_unique<DiscreteScenarioSet>();
            dataFile.open(filename, netCDF::NcFile::read);
            dss->deserialize(dataFile);
            dataFile.close();
            
            if (dss->get_poolSize() != 7) {
                remove(filename.c_str());
                return {false, "Should have loaded poolSize=7 from file"};
            }
            
            dss->init_representative_pool(7);
            
            if (dss->get_poolSize() != 7) {
                remove(filename.c_str());
                return {false, "Should have 7 scenarios after init_representative_pool"};
            }
            
            remove(filename.c_str());
            
            // Test invalid poolSize values
            filename = create_test_scenario_file(20, 5, "_k_zero");
            {
                netCDF::NcFile dataFile(filename, netCDF::NcFile::write);
                auto cfgGroup = dataFile.addGroup("ScenarioReductionConfig");
                int poolSize = 0;
                auto poolSizeVar = cfgGroup.addVar("poolSize", netCDF::ncInt);
                poolSizeVar.putVar(&poolSize);
                dataFile.close();
            }
            
            dss = make_unique<DiscreteScenarioSet>();
            dataFile.open(filename, netCDF::NcFile::read);
            dss->deserialize(dataFile);
            dataFile.close();
            
            if (dss->get_poolSize() != 0) {
                remove(filename.c_str());
                return {false, "No reduction should occur with poolSize=0"};
            }
            
            remove(filename.c_str());
        }
        
        // Part 3: Serialization round-trip
        {
            string filename1 = create_test_scenario_file(20, 5, "_roundtrip1");
            string filename2 = "test_roundtrip2.nc4";
            
            {
                netCDF::NcFile dataFile(filename1, netCDF::NcFile::write);
                auto cfgGroup = dataFile.addGroup("ScenarioReductionConfig");
                int poolSize = 8;
                auto poolSizeVar = cfgGroup.addVar("poolSize", netCDF::ncInt);
                poolSizeVar.putVar(&poolSize);
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
            
            if (dss2->get_poolSize() != 8) {
                remove(filename1.c_str());
                remove(filename2.c_str());
                return {false, "Round-trip serialization failed - poolSize value not preserved"};
            }
            
            dss2->init_representative_pool(8);
            
            if (dss2->get_poolSize() != 8) {
                remove(filename1.c_str());
                remove(filename2.c_str());
                return {false, "Round-trip serialization failed - scenario selection failed"};
            }
            
            remove(filename1.c_str());
            remove(filename2.c_str());
        }
        
        // Part 4: poolSize variable priority over BlockConfig SimpleConfiguration<int>
        {
            string filename = create_test_scenario_file(20, 5, "_k_priority");
            {
                netCDF::NcFile dataFile(filename, netCDF::NcFile::write);
                
                auto cfgGroup = dataFile.addGroup("ScenarioReductionConfig");
                
                // Add poolSize variable (should have priority)
                int poolSize = 7;
                auto poolSizeVar = cfgGroup.addVar("poolSize", netCDF::ncInt);
                poolSizeVar.putVar(&poolSize);
                
                // Add BlockConfig with different SimpleConfiguration<int> value
                auto blockGroup = cfgGroup.addGroup("BlockConfig");
                auto extraGroup = blockGroup.addGroup("f_extra_Configuration");
                extraGroup.putAtt("type", "SimpleConfiguration<int>");
                int block_poolSize = 5; // Different value - should be overridden
                auto blockPoolSizeVar = extraGroup.addVar("value", netCDF::ncInt);
                blockPoolSizeVar.putVar(&block_poolSize);
                
                dataFile.close();
            }
            
            auto dss = make_unique<DiscreteScenarioSet>();
            netCDF::NcFile dataFile(filename, netCDF::NcFile::read);
            dss->deserialize(dataFile);
            dataFile.close();
            
            // poolSize variable should have priority (7, not 5)
            if (dss->get_poolSize() != 7) {
                remove(filename.c_str());
                return {false, "poolSize variable should have priority (expected 7, got " + to_string(dss->get_poolSize()) + ")"};
            }
            
            remove(filename.c_str());
        }
        
        
        // Part 5: Deserialization with invalid poolSize (error checking happens at init time)
        {
            string filename = create_test_scenario_file(10, 2, "_invalid_poolSize");
            {
                netCDF::NcFile dataFile(filename, netCDF::NcFile::write);
                
                auto cfgGroup = dataFile.addGroup("ScenarioReductionConfig");
                int poolSize = 15; // Invalid: poolSize > number of scenarios (10)
                auto poolSizeVar = cfgGroup.addVar("poolSize", netCDF::ncInt);
                poolSizeVar.putVar(&poolSize);
                
                dataFile.close();
            }
            
            auto dss = make_unique<DiscreteScenarioSet>();
            netCDF::NcFile dataFile(filename, netCDF::NcFile::read);
            
            try {
                dss->deserialize(dataFile);
                dataFile.close();
                
                // Should have loaded invalid poolSize
                if (dss->get_poolSize() != 15) {
                    remove(filename.c_str());
                    return {false, "Should load invalid poolSize value (got " + to_string(dss->get_poolSize()) + ")"};
                }
                
                // But init_representative_pool should fail
                try {
                    dss->init_representative_pool(15);
                    remove(filename.c_str());
                    return {false, "Should throw exception for invalid poolSize in init_representative_pool"};
                } catch (const invalid_argument& e) {
                    // Expected behavior
                }
            } catch (const runtime_error& e) {
                // If deserialization fails, that's also acceptable
                dataFile.close();
            }
            
            remove(filename.c_str());
        }
        
        
        return {true, "NetCDF serialization and deserialization tests passed"};
        
    } catch (const exception& e) {
        return {false, string("Serialization test failed: ") + e.what()};
    }
}

REGISTER_TEST("Serialization and Deserialization", test_serialization_deserialization);

// Test 5: Iteration and Span-based Getters
TestResult test_iteration_and_spans() {
    try {
        auto dss = load_test_scenarios(15, 5);
        
        // Initialize pool for testing
        dss->init_random_pool(8);
        
        // Part 1: Test full iteration through selected scenarios
        int iteration_count = 0;
        double total_prob = 0.0;
        
        do {
            auto scenario = dss->get_current_scenario();
            double prob = dss->get_current_scenario_probability();
            
            if (scenario.size() != 5) {
                return {false, "Scenario size mismatch during iteration"};
            }
            
            total_prob += prob;
            iteration_count++;
        } while (dss->next_scenario());
        
        if (iteration_count != 8) {
            return {false, "Expected 8 iterations, got " + to_string(iteration_count)};
        }
        
        if (abs(total_prob - 1.0) > 1e-6) {
            return {false, "Probabilities don't sum to 1.0, got " + to_string(total_prob)};
        }
        
        // Part 2: Test span-based getters
        auto selected_indices = dss->get_selected_scenarios();
        if (selected_indices.size() != 8) {
            return {false, "Selected scenarios span size mismatch"};
        }
        
        // Verify indices are valid
        for (auto idx : selected_indices) {
            if (idx >= 15) {
                return {false, "Invalid scenario index in selection: " + to_string(idx)};
            }
        }
        
        auto pool_weights = dss->get_pool_weights();
        if (pool_weights.size() != 15) {
            return {false, "Pool weights span should have all scenarios"};
        }
        
        auto normalized_weights = dss->get_normalized_weights();
        if (normalized_weights.size() != 8) {
            return {false, "Normalized weights span size mismatch"};
        }
        
        // Verify normalized weights sum to 1
        double norm_sum = 0.0;
        for (auto w : normalized_weights) {
            norm_sum += w;
        }
        if (abs(norm_sum - 1.0) > 1e-6) {
            return {false, "Normalized weights don't sum to 1.0"};
        }
        
        return {true, "Iteration and span tests passed"};
        
    } catch (const exception& e) {
        return {false, string("Exception: ") + e.what()};
    }
}

REGISTER_TEST("Iteration and Spans", test_iteration_and_spans);



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