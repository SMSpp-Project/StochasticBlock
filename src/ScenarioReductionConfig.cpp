/*--------------------------------------------------------------------------*/
/*------------------- File ScenarioReductionConfig.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the ScenarioReductionConfig class.
 *
 * \author Benoît Tran \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Benoît Tran
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ScenarioReductionConfig.h"
#include "SimpleConfiguration.h"
#include "ComputeConfig.h"
#include "MILPSolver.h"


/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------- FACTORY MANAGEMENT ----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0(ScenarioReductionConfig);

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

void ScenarioReductionConfig::initialize(int k, float ell, const std::string& algorithm)
{
    // Validate parameters
    validate_k(k);
    validate_ell(ell);
    validate_algorithm(algorithm);
    
    // Get or create k, ell and algorithm members and set their values
    get_or_create_cfl_member("k", k)->set_value(k);
    get_or_create_cfl_member("ell", ell)->set_value(ell);
    get_or_create_solver_member("algorithm", algorithm)->set_value(algorithm);
}

ScenarioReductionConfig::ScenarioReductionConfig(int k)
    : Configuration(), f_cfl_config(nullptr), f_solver_config(nullptr)
{
    // Initialize with k and default values for ell and algorithm
    initialize(k, DEFAULT_ELL, DEFAULT_ALGORITHM);
}

/*--------------------------------------------------------------------------*/

ScenarioReductionConfig::ScenarioReductionConfig(const netCDF::NcGroup& group)
    : Configuration(), f_cfl_config(nullptr), f_solver_config(nullptr)
{
    // Directly deserialize from the netCDF group
    // This will create and populate the config objects as needed
    ScenarioReductionConfig::deserialize(group);
}

/*--------------------------------------------------------------------------*/

ScenarioReductionConfig::ScenarioReductionConfig(int k, float ell, const std::string& algorithm)
    : Configuration(), f_cfl_config(nullptr), f_solver_config(nullptr)
{
    // Initialize with provided parameters
    initialize(k, ell, algorithm);
}

/*--------------------------------------------------------------------------*/

ScenarioReductionConfig* ScenarioReductionConfig::clone() const
{
    // Create a new empty ScenarioReductionConfig
    auto* config = new ScenarioReductionConfig(DEFAULT_K);
    
    try {
        // Clone the CFLConfig if it exists
        if (f_cfl_config) {
            delete config->f_cfl_config;
            config->f_cfl_config = static_cast<BlockConfig*>(f_cfl_config->clone());
        } else {
            // Create a new CFLConfig with default values if the original is null
            config->set_k(DEFAULT_K);
            config->set_ell(DEFAULT_ELL);
        }
        
        // Clone the SolverConfig if it exists
        if (f_solver_config) {
            delete config->f_solver_config;
            config->f_solver_config = static_cast<BlockSolverConfig*>(f_solver_config->clone());
        } else {
            // Create a new SolverConfig with default algorithm if the original is null
            config->set_algorithm(DEFAULT_ALGORITHM);
        }
        
        // Validate the cloned configuration
        if (config->f_cfl_config) {
            validate_cfl_config(config->f_cfl_config, false);  // Don't require all members in validation
        }
        
        if (config->f_solver_config) {
            validate_solver_config(config->f_solver_config, false);  // Don't require all members in validation
        }
        
    } catch (const std::exception& e) {
        // If anything fails during cloning, clean up and rethrow
        delete config;
        throw std::runtime_error("Error cloning ScenarioReductionConfig: " + std::string(e.what()));
    }
    
    return config;
}

/*--------------------------------------------------------------------------*/

void ScenarioReductionConfig::reset()
{
    // Delete existing configurations
    delete f_cfl_config;
    delete f_solver_config;
    
    // Create new configurations with default values
    f_cfl_config = new BlockConfig();
    f_cfl_config->add_member("k", new SimpleConfiguration<int>(DEFAULT_K));
    f_cfl_config->add_member("ell", new SimpleConfiguration<float>(DEFAULT_ELL));
    
    f_solver_config = new BlockSolverConfig(false);
    f_solver_config->add_member("algorithm", new SimpleConfiguration<std::string>(DEFAULT_ALGORITHM));
    
    report_message("Configuration reset to default values", 0);
}

/*--------------------------------------------------------------------------*/

bool ScenarioReductionConfig::check_consistency(bool auto_fix)
{
    bool is_consistent = true;
    std::vector<std::string> issues;
    
    // Check CFLConfig
    if (!f_cfl_config) {
        issues.push_back("CFLConfig is null");
        is_consistent = false;
        
        if (auto_fix) {
            report_message("Creating missing CFLConfig", 0);
            f_cfl_config = new BlockConfig();
        }
    }
    
    if (f_cfl_config) {
        // Check k parameter
        if (!f_cfl_config->has_member("k")) {
            issues.push_back("Missing k parameter in CFLConfig");
            is_consistent = false;
            
            if (auto_fix) {
                report_message("Adding missing k parameter with default value", 0);
                f_cfl_config->add_member("k", new SimpleConfiguration<int>(DEFAULT_K));
            }
        } else {
            auto* k_config = dynamic_cast<SimpleConfiguration<int>*>(f_cfl_config->get_member("k"));
            if (!k_config) {
                issues.push_back("k parameter has wrong type in CFLConfig");
                is_consistent = false;
                
                if (auto_fix) {
                    report_message("Replacing k parameter with default value", 0);
                    delete f_cfl_config->get_member("k");
                    f_cfl_config->add_member("k", new SimpleConfiguration<int>(DEFAULT_K));
                }
            } else {
                // Check k parameter value
                int k = k_config->get_value();
                if (k <= 0) {
                    issues.push_back("k parameter is not positive");
                    is_consistent = false;
                    
                    if (auto_fix) {
                        report_message("Setting k parameter to default value", 0);
                        k_config->set_value(DEFAULT_K);
                    }
                }
            }
        }
        
        // Check ell parameter
        if (!f_cfl_config->has_member("ell")) {
            issues.push_back("Missing ell parameter in CFLConfig");
            is_consistent = false;
            
            if (auto_fix) {
                report_message("Adding missing ell parameter with default value", 0);
                f_cfl_config->add_member("ell", new SimpleConfiguration<float>(DEFAULT_ELL));
            }
        } else {
            auto* ell_config = dynamic_cast<SimpleConfiguration<float>*>(f_cfl_config->get_member("ell"));
            if (!ell_config) {
                issues.push_back("ell parameter has wrong type in CFLConfig");
                is_consistent = false;
                
                if (auto_fix) {
                    report_message("Replacing ell parameter with default value", 0);
                    delete f_cfl_config->get_member("ell");
                    f_cfl_config->add_member("ell", new SimpleConfiguration<float>(DEFAULT_ELL));
                }
            } else {
                // Check ell parameter value
                float ell = ell_config->get_value();
                if (ell <= 0.0f) {
                    issues.push_back("ell parameter is not positive");
                    is_consistent = false;
                    
                    if (auto_fix) {
                        report_message("Setting ell parameter to default value", 0);
                        ell_config->set_value(DEFAULT_ELL);
                    }
                }
            }
        }
    }
    
    // Check SolverConfig
    if (!f_solver_config) {
        issues.push_back("SolverConfig is null");
        is_consistent = false;
        
        if (auto_fix) {
            report_message("Creating missing SolverConfig", 0);
            f_solver_config = new BlockSolverConfig(false);
        }
    }
    
    if (f_solver_config) {
        // Check algorithm parameter
        if (!f_solver_config->has_member("algorithm")) {
            issues.push_back("Missing algorithm parameter in SolverConfig");
            is_consistent = false;
            
            if (auto_fix) {
                report_message("Adding missing algorithm parameter with default value", 0);
                f_solver_config->add_member("algorithm", new SimpleConfiguration<std::string>(DEFAULT_ALGORITHM));
            }
        } else {
            auto* algo_config = dynamic_cast<SimpleConfiguration<std::string>*>(f_solver_config->get_member("algorithm"));
            if (!algo_config) {
                issues.push_back("algorithm parameter has wrong type in SolverConfig");
                is_consistent = false;
                
                if (auto_fix) {
                    report_message("Replacing algorithm parameter with default value", 0);
                    delete f_solver_config->get_member("algorithm");
                    f_solver_config->add_member("algorithm", new SimpleConfiguration<std::string>(DEFAULT_ALGORITHM));
                }
            } else {
                // Check algorithm parameter value
                std::string algorithm = algo_config->get_value();
                if (!is_valid_algorithm(algorithm)) {
                    issues.push_back("Invalid algorithm: " + algorithm);
                    is_consistent = false;
                    
                    if (auto_fix) {
                        report_message("Setting algorithm to default value", 0);
                        algo_config->set_value(DEFAULT_ALGORITHM);
                    }
                }
            }
        }
        
        // Check compute config consistency with algorithm
        if (f_solver_config->has_ComputeConfig(0)) {
            std::string computeName = f_solver_config->get_ComputeName(0);
            std::string algorithm = get_algorithm();
            
            if (algorithm == "MILP" && computeName != "HiGHSMILPSolver") {
                issues.push_back("Algorithm 'MILP' requires HiGHSMILPSolver, but " + computeName + " is set");
                is_consistent = false;
                
                if (auto_fix) {
                    report_message("Setting ComputeName to HiGHSMILPSolver for MILP algorithm", 0);
                    f_solver_config->set_ComputeName(0, "HiGHSMILPSolver");
                }
            } else if (algorithm != "MILP" && computeName != "ScenarioReductionSolver") {
                issues.push_back("Non-MILP algorithms require ScenarioReductionSolver, but " + computeName + " is set");
                is_consistent = false;
                
                if (auto_fix) {
                    report_message("Setting ComputeName to ScenarioReductionSolver for " + algorithm + " algorithm", 0);
                    f_solver_config->set_ComputeName(0, "ScenarioReductionSolver");
                }
            }
        }
    }
    
    // Report all issues
    if (!issues.empty()) {
        for (const auto& issue : issues) {
            report_message("Consistency issue: " + issue, auto_fix ? 0 : 2);
        }
    }
    
    return is_consistent;
}

/*--------------------------------------------------------------------------*/

ScenarioReductionConfig::~ScenarioReductionConfig()
{
    delete f_cfl_config;
    delete f_solver_config;
}

/*--------------------------------------------------------------------------*/
/*------------------------- SERIALIZATION METHODS -------------------------*/
/*--------------------------------------------------------------------------*/

BlockConfig* ScenarioReductionConfig::deserialize_cfl_config(const netCDF::NcGroup& group) const
{
    // Look for CFLConfig subgroup - this is mandatory
    netCDF::NcGroup cflConfig = group.getGroup("CFLConfig");
    if (cflConfig.isNull()) {
        throw std::invalid_argument("ScenarioReductionConfig: Missing mandatory CFLConfig group");
    }
    
    // Read k parameter (mandatory)
    int k;
    netCDF::NcGroupAtt kAtt = cflConfig.getAtt("k");
    if (!kAtt.isNull()) {
        kAtt.getValues(&k);
        validate_k(k);
    } else {
        throw std::invalid_argument("ScenarioReductionConfig: Missing mandatory 'k' attribute in CFLConfig");
    }
    
    // Read ell parameter (optional with default 2.0)
    float ell = DEFAULT_ELL;
    netCDF::NcGroupAtt ellAtt = cflConfig.getAtt("ell");
    if (!ellAtt.isNull()) {
        ellAtt.getValues(&ell);
        validate_ell(ell);
    }
    
    // Create and return the new BlockConfig
    auto* config = new BlockConfig();
    config->add_member("k", new SimpleConfiguration<int>(k));
    config->add_member("ell", new SimpleConfiguration<float>(ell));
    
    return config;
}

/*--------------------------------------------------------------------------*/

BlockSolverConfig* ScenarioReductionConfig::deserialize_solver_config(const netCDF::NcGroup& group) const
{
    // Look for SolverConfig subgroup - this is mandatory
    netCDF::NcGroup solverConfig = group.getGroup("SolverConfig");
    if (solverConfig.isNull()) {
        throw std::invalid_argument("ScenarioReductionConfig: Missing mandatory SolverConfig group");
    }
    
    // Check if it contains a complete BlockSolverConfig
    netCDF::NcGroupAtt typeAtt = solverConfig.getAtt("type");
    if (!typeAtt.isNull()) {
        std::string type;
        typeAtt.getValues(type);
        
        // If it's a proper BlockSolverConfig, deserialize it
        if (type == "BlockSolverConfig" || type == "RBlockSolverConfig") {
            auto* config = dynamic_cast<BlockSolverConfig*>(Configuration::new_Configuration(solverConfig));
            
            // Check for algorithm
            if (!config->has_member("algorithm")) {
                throw std::invalid_argument("ScenarioReductionConfig: Missing 'algorithm' in SolverConfig");
            }
            
            return config;
        }
    }
    
    // Otherwise, read the algorithm directly from an attribute
    std::string algorithm;
    netCDF::NcGroupAtt algoAtt = solverConfig.getAtt("algorithm");
    if (!algoAtt.isNull()) {
        algoAtt.getValues(algorithm);
    } else {
        throw std::invalid_argument("ScenarioReductionConfig: Missing mandatory 'algorithm' attribute in SolverConfig");
    }
    
    // Validate the algorithm
    validate_algorithm(algorithm);
    
    // Create and return a new SolverConfig with the algorithm
    auto* config = new BlockSolverConfig(false);
    config->add_member("algorithm", new SimpleConfiguration<std::string>(algorithm));
    
    return config;
}

/*--------------------------------------------------------------------------*/

void ScenarioReductionConfig::deserialize(const netCDF::NcGroup& group)
{
    // Call the base class implementation first
    Configuration::deserialize(group);
    
    // Create placeholder configs if they don't exist
    ensure_cfl_config();
    ensure_solver_config();
    
    try {
        // Process CFLConfig
        BlockConfig* new_cfl_config = deserialize_cfl_config(group);
        delete f_cfl_config;
        f_cfl_config = new_cfl_config;
        
        // Process SolverConfig
        BlockSolverConfig* new_solver_config = deserialize_solver_config(group);
        delete f_solver_config;
        f_solver_config = new_solver_config;
        
    } catch (const std::exception& e) {
        throw std::invalid_argument("ScenarioReductionConfig: Error during deserialization: " + 
                                   std::string(e.what()));
    }
}

/*--------------------------------------------------------------------------*/

void ScenarioReductionConfig::serialize(netCDF::NcGroup& group) const
{
    // Call the base class implementation first
    Configuration::serialize(group);
    
    // Create CFLConfig subgroup
    netCDF::NcGroup cflConfig = group.addGroup("CFLConfig");
    
    // Get k from CFLConfig
    int k = get_k();
    cflConfig.putAtt("k", netCDF::NcInt(), k);
    
    // Get ell from CFLConfig
    float ell = get_ell();
    cflConfig.putAtt("ell", netCDF::NcFloat(), ell);
    
    // Serialize the CFLConfig if it exists
    if (f_cfl_config) {
        f_cfl_config->serialize(cflConfig);
    }
    
    // Create SolverConfig subgroup
    netCDF::NcGroup solverConfig = group.addGroup("SolverConfig");
    
    // Get algorithm from SolverConfig
    std::string algorithm = get_algorithm();
    solverConfig.putAtt("algorithm", algorithm);
    
    // Serialize the SolverConfig if it exists
    if (f_solver_config) {
        f_solver_config->serialize(solverConfig);
    }
}

/*--------------------------------------------------------------------------*/
/*-------------------------- ACCESSORS/MUTATORS ---------------------------*/
/*--------------------------------------------------------------------------*/

void ScenarioReductionConfig::set_k(int value)
{
    // Validate the parameter
    validate_k(value);
    
    // Get or create k member and update its value
    get_or_create_cfl_member("k", value)->set_value(value);
}

/*--------------------------------------------------------------------------*/

int ScenarioReductionConfig::get_k() const
{
    // Check if CFLConfig exists
    if (!f_cfl_config) {
        report_using_default("k", std::to_string(DEFAULT_K));
        return DEFAULT_K;
    }
    
    // Check if k member exists
    if (!f_cfl_config->has_member("k")) {
        report_using_default("k", std::to_string(DEFAULT_K));
        return DEFAULT_K;
    }
    
    // Check if k member has correct type
    auto* k_config = dynamic_cast<SimpleConfiguration<int>*>(f_cfl_config->get_member("k"));
    if (!k_config) {
        report_using_default("k", std::to_string(DEFAULT_K));
        return DEFAULT_K;
    }
    
    return k_config->get_value();
}

/*--------------------------------------------------------------------------*/

void ScenarioReductionConfig::set_ell(float value)
{
    // Validate the parameter
    validate_ell(value);
    
    // Get or create ell member and update its value
    get_or_create_cfl_member("ell", value)->set_value(value);
}

/*--------------------------------------------------------------------------*/

float ScenarioReductionConfig::get_ell() const
{
    // Check if CFLConfig exists
    if (!f_cfl_config) {
        report_using_default("ell", std::to_string(DEFAULT_ELL));
        return DEFAULT_ELL;
    }
    
    // Check if ell member exists
    if (!f_cfl_config->has_member("ell")) {
        report_using_default("ell", std::to_string(DEFAULT_ELL));
        return DEFAULT_ELL;
    }
    
    // Check if ell member has correct type
    auto* ell_config = dynamic_cast<SimpleConfiguration<float>*>(f_cfl_config->get_member("ell"));
    if (!ell_config) {
        report_using_default("ell", std::to_string(DEFAULT_ELL));
        return DEFAULT_ELL;
    }
    
    return ell_config->get_value();
}

/*--------------------------------------------------------------------------*/

void ScenarioReductionConfig::set_algorithm(const std::string& alg)
{
    // Validate the algorithm choice
    validate_algorithm(alg);
    
    // Get or create algorithm member and update its value
    get_or_create_solver_member("algorithm", alg)->set_value(alg);
}

/*--------------------------------------------------------------------------*/

const std::string& ScenarioReductionConfig::get_algorithm() const
{
    // Check if SolverConfig exists
    if (!f_solver_config) {
        report_using_default("algorithm", DEFAULT_ALGORITHM);
        return DEFAULT_ALGORITHM;
    }
    
    // Check if algorithm member exists
    if (!f_solver_config->has_member("algorithm")) {
        report_using_default("algorithm", DEFAULT_ALGORITHM);
        return DEFAULT_ALGORITHM;
    }
    
    // Check if algorithm member has correct type
    auto* algo_config = dynamic_cast<SimpleConfiguration<std::string>*>(f_solver_config->get_member("algorithm"));
    if (!algo_config) {
        report_using_default("algorithm", DEFAULT_ALGORITHM);
        return DEFAULT_ALGORITHM;
    }
    
    return algo_config->get_value();
}

/*--------------------------------------------------------------------------*/

void ScenarioReductionConfig::set_solver_config(BlockSolverConfig* config)
{
    // Check if we need to add a default algorithm member
    if (config && !config->has_member("algorithm")) {
        report_message("Adding missing algorithm member to supplied SolverConfig");
        config->add_member("algorithm", new SimpleConfiguration<std::string>(DEFAULT_ALGORITHM));
    }
    
    // Validate the configuration (will throw an exception if invalid)
    if (config) {
        validate_solver_config(config);
    }
    
    // Replace the old configuration with the new one
    delete f_solver_config;
    f_solver_config = config;
}

/*--------------------------------------------------------------------------*/

BlockSolverConfig* ScenarioReductionConfig::get_solver_config() const
{
    if (!f_solver_config) {
        report_message("SolverConfig is null in get_solver_config");
    }
    return f_solver_config;
}

/*--------------------------------------------------------------------------*/
/*----------------------- CONFIGURATION EXTRACTION ------------------------*/
/*--------------------------------------------------------------------------*/

BlockConfig* ScenarioReductionConfig::extract_block_config(size_t n_scenarios) const
{
    // If we have a valid CFLConfig, clone it, otherwise create a new one
    BlockConfig* blockConfig;
    if (f_cfl_config) {
        blockConfig = static_cast<BlockConfig*>(f_cfl_config->clone());
    } else {
        report_message("Creating new BlockConfig in extract_block_config because CFLConfig is null");
        blockConfig = new BlockConfig();
        
        // Set defaults
        blockConfig->add_member("k", new SimpleConfiguration<int>(DEFAULT_K));
        blockConfig->add_member("ell", new SimpleConfiguration<float>(DEFAULT_ELL));
    }
    
    // Make sure the required parameters are present with correct values
    int k = get_k();
    float ell = get_ell();
    
    // Helper function to update or add a configuration member
    auto update_or_add_member = [this, &blockConfig](const std::string& name, auto value) {
        using T = decltype(value);
        if (blockConfig->has_member(name)) {
            auto* config = dynamic_cast<SimpleConfiguration<T>*>(blockConfig->get_member(name));
            if (config) {
                config->set_value(value);
            } else {
                report_message("Replacing incompatible " + name + " member in extracted BlockConfig");
                delete blockConfig->get_member(name);
                blockConfig->add_member(name, new SimpleConfiguration<T>(value));
            }
        } else {
            blockConfig->add_member(name, new SimpleConfiguration<T>(value));
        }
    };
    
    // Update or add k, ell, and n parameters
    update_or_add_member("k", k);
    update_or_add_member("ell", ell);
    update_or_add_member("n", static_cast<int>(n_scenarios));
    
    return blockConfig;
}

/*--------------------------------------------------------------------------*/

BlockSolverConfig* ScenarioReductionConfig::extract_solver_config() const
{
    // Get the algorithm we're using
    std::string algorithm = get_algorithm();
    float ell = get_ell();
    
    // If we have a valid SolverConfig, clone it
    if (f_solver_config) {
        auto* solverConfig = static_cast<BlockSolverConfig*>(f_solver_config->clone());
        
        // If it doesn't have all the necessary configurations, add them
        if (!solverConfig->has_ComputeConfig(0)) {
            if (algorithm == "MILP") {
                // Add MILPSolver configuration
                solverConfig->add_ComputeConfig("HiGHSMILPSolver");
                
                // Create a ComputeConfig for the MILPSolver
                auto* computeConfig = new ComputeConfig();
                computeConfig->set_int_par("intLogVerb", 0);  // Low verbosity
                computeConfig->set_dbl_par("dblRelAcc", 1e-6);  // Relative accuracy
                solverConfig->set_ComputeConfig(0, computeConfig);
                
                report_message("Added missing MILPSolver configuration to cloned SolverConfig");
            } else {
                // Add ScenarioReductionSolver configuration
                solverConfig->add_ComputeConfig("ScenarioReductionSolver");
                
                // Create a ComputeConfig for the solver
                auto* computeConfig = new ComputeConfig();
                computeConfig->set_str_par("algorithm", algorithm);
                computeConfig->set_dbl_par("ell", static_cast<double>(ell));
                solverConfig->set_ComputeConfig(0, computeConfig);
                
                report_message("Added missing ScenarioReductionSolver configuration to cloned SolverConfig");
            }
        } else {
            // Make sure the compute config has the correct parameters
            auto* computeConfig = solverConfig->get_ComputeConfig(0);
            if (computeConfig) {
                if (algorithm == "MILP") {
                    // Check if this is the correct solver type
                    if (solverConfig->get_ComputeName(0) != "HiGHSMILPSolver") {
                        report_message("Incorrect solver type in SolverConfig, replacing with HiGHSMILPSolver");
                        solverConfig->set_ComputeName(0, "HiGHSMILPSolver");
                    }
                    
                    // Set default parameters if needed
                    if (!computeConfig->has_int_par("intLogVerb")) {
                        computeConfig->set_int_par("intLogVerb", 0);
                    }
                    if (!computeConfig->has_dbl_par("dblRelAcc")) {
                        computeConfig->set_dbl_par("dblRelAcc", 1e-6);
                    }
                } else {
                    // Check if this is the correct solver type
                    if (solverConfig->get_ComputeName(0) != "ScenarioReductionSolver") {
                        report_message("Incorrect solver type in SolverConfig, replacing with ScenarioReductionSolver");
                        solverConfig->set_ComputeName(0, "ScenarioReductionSolver");
                    }
                    
                    // Update algorithm
                    computeConfig->set_str_par("algorithm", algorithm);
                    
                    // Update ell parameter
                    computeConfig->set_dbl_par("ell", static_cast<double>(ell));
                }
            }
        }
        
        return solverConfig;
    }

    // If we don't have a SolverConfig, create a new one
    report_message("Creating new BlockSolverConfig in extract_solver_config because SolverConfig is null");
    auto* solverConfig = new BlockSolverConfig(false);  // Not differential

    // Determine which solver to use based on the algorithm
    if (algorithm == "MILP") {
        // Add a MILPSolver
        solverConfig->add_ComputeConfig("HiGHSMILPSolver");
        
        // Create a ComputeConfig for the MILPSolver
        auto* computeConfig = new ComputeConfig();
        
        // Set solver parameters
        computeConfig->set_int_par("intLogVerb", 0);  // Low verbosity
        computeConfig->set_dbl_par("dblRelAcc", 1e-6);  // Relative accuracy
        
        // Add the ComputeConfig to the BlockSolverConfig
        solverConfig->set_ComputeConfig(0, computeConfig);
    } else {
        // Use ScenarioReductionSolver with the specified algorithm
        solverConfig->add_ComputeConfig("ScenarioReductionSolver");
        
        // Create a ComputeConfig for the solver
        auto* computeConfig = new ComputeConfig();
        
        // Set the algorithm
        computeConfig->set_str_par("algorithm", algorithm);
        
        // Set the ell parameter
        computeConfig->set_dbl_par("ell", static_cast<double>(ell));
        
        // Add the ComputeConfig to the BlockSolverConfig
        solverConfig->set_ComputeConfig(0, computeConfig);
    }

    return solverConfig;
}

/*--------------------------------------------------------------------------*/
/*-------------------------- UTILITY METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

const std::vector<std::string>& ScenarioReductionConfig::get_valid_algorithms()
{
    // This static vector contains string representations of ScenarioReductionSolver::Algorithm enum values
    static const std::vector<std::string> valid_algorithms = {
        "Baseline",  // ScenarioReductionSolver::Algorithm::Baseline
        "Dupacova",  // ScenarioReductionSolver::Algorithm::Dupacova  
        "BestFit",   // ScenarioReductionSolver::Algorithm::BestFit
        "FirstFit",  // ScenarioReductionSolver::Algorithm::FirstFit
        "MILP"       // ScenarioReductionSolver::Algorithm::MILP
    };
    
    return valid_algorithms;
}

/*--------------------------------------------------------------------------*/

bool ScenarioReductionConfig::is_valid_algorithm(const std::string& algorithm)
{
    const auto& valid_algorithms = get_valid_algorithms();
    return std::find(valid_algorithms.begin(), valid_algorithms.end(), algorithm) != valid_algorithms.end();
}

/*--------------------------------------------------------------------------*/

std::string ScenarioReductionConfig::get_valid_algorithm_options()
{
    const auto& valid_algorithms = get_valid_algorithms();
    if (valid_algorithms.empty()) {
        return "";
    }
    
    std::string options = valid_algorithms[0];
    for (size_t i = 1; i < valid_algorithms.size(); ++i) {
        options += ", " + valid_algorithms[i];
    }
    
    return options;
}

/*--------------------------------------------------------------------------*/
/*----------------------------- I/O METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void ScenarioReductionConfig::print(std::ostream& output) const
{
    output << "ScenarioReductionConfig [" << this << "]" << std::endl;
    
    // Print CFLConfig details
    output << "  CFLConfig:" << std::endl;
    if (f_cfl_config) {
        output << "    k: " << get_k() << std::endl;
        output << "    ell: " << get_ell() << std::endl;
    } else {
        output << "    <null>" << std::endl;
    }
    
    // Print SolverConfig details
    output << "  SolverConfig:" << std::endl;
    if (f_solver_config) {
        output << "    algorithm: " << get_algorithm() << std::endl;
        output << "    type: " << f_solver_config->classname() << " [" << f_solver_config << "]" << std::endl;
        
        // Print ComputeConfig details if available
        if (f_solver_config->has_ComputeConfig(0)) {
            output << "    ComputeConfig[0]: " << f_solver_config->get_ComputeName(0) << std::endl;
        }
    } else {
        output << "    <null>" << std::endl;
    }
}

/*--------------------------------------------------------------------------*/

void ScenarioReductionConfig::load(std::istream& input)
{
    // Create empty configs if they don't exist
    if (!f_cfl_config) {
        f_cfl_config = new BlockConfig();
    }
    
    if (!f_solver_config) {
        f_solver_config = new BlockSolverConfig(false);
    }
    
    // Read parameters for CFLConfig
    int k;
    float ell;
    input >> eatcomments >> k;
    input >> eatcomments >> ell;
    
    // Validate parameters
    validate_k(k);
    validate_ell(ell);
    
    // Update CFLConfig
    set_k(k);
    set_ell(ell);
    
    // Read algorithm for SolverConfig
    std::string algorithm;
    input >> eatcomments >> algorithm;
    set_algorithm(algorithm);
    
    // Check if there's a BlockSolverConfig to load
    std::string configType;
    input >> eatcomments >> configType;
    
    if (configType != "*") {
        // A configuration type is specified, create and load it
        delete f_solver_config;
        f_solver_config = dynamic_cast<BlockSolverConfig*>(Configuration::new_Configuration(configType));
        if (f_solver_config) {
            input >> *f_solver_config;
            
            // Make sure it has an algorithm member
            if (!f_solver_config->has_member("algorithm")) {
                std::cerr << "Warning: Adding missing algorithm member to loaded SolverConfig" << std::endl;
                f_solver_config->add_member("algorithm", new SimpleConfiguration<std::string>(algorithm));
            }
        }
    }
}

/*--------------------------------------------------------------------------*/
/*--------------------- HELPER METHODS FOR CONFIG MEMBERS --------------------*/
/*--------------------------------------------------------------------------*/

BlockConfig* ScenarioReductionConfig::ensure_cfl_config()
{
    if (!f_cfl_config) {
        report_message("Creating a new CFLConfig");
        f_cfl_config = new BlockConfig();
    }
    return f_cfl_config;
}

/*--------------------------------------------------------------------------*/

BlockSolverConfig* ScenarioReductionConfig::ensure_solver_config()
{
    if (!f_solver_config) {
        report_message("Creating a new SolverConfig");
        f_solver_config = new BlockSolverConfig(false);
    }
    return f_solver_config;
}

/*--------------------------------------------------------------------------*/

void ScenarioReductionConfig::report_message(const std::string& message, int severity)
{
    // Determine the prefix based on severity
    std::string prefix;
    switch (severity) {
        case 0:
            prefix = "Info: ";
            break;
        case 1:
            prefix = "Warning: ";
            break;
        case 2:
            prefix = "Error: ";
            break;
        default:
            prefix = "Note: ";
    }
    
    // Output to stderr with class name and message
    std::cerr << "ScenarioReductionConfig: " << prefix << message << std::endl;
}

/*--------------------------------------------------------------------------*/

void ScenarioReductionConfig::report_using_default(const std::string& param_name, const std::string& default_value)
{
    report_message(param_name + " missing or invalid, using default value: " + default_value);
}

/*--------------------------------------------------------------------------*/

void ScenarioReductionConfig::validate_k(int k)
{
    if (k <= 0) {
        throw std::invalid_argument("ScenarioReductionConfig: k must be positive");
    }
}

/*--------------------------------------------------------------------------*/

void ScenarioReductionConfig::validate_ell(float ell)
{
    if (ell <= 0.0f) {
        throw std::invalid_argument("ScenarioReductionConfig: ell must be positive");
    }
}

/*--------------------------------------------------------------------------*/

void ScenarioReductionConfig::validate_algorithm(const std::string& algorithm)
{
    if (!is_valid_algorithm(algorithm)) {
        throw std::invalid_argument("ScenarioReductionConfig: Invalid algorithm '" + 
                                   algorithm + "'. Must be one of: " + get_valid_algorithm_options());
    }
}

/*--------------------------------------------------------------------------*/

void ScenarioReductionConfig::validate_cfl_config(const BlockConfig* config, bool check_required)
{
    if (!config) {
        throw std::invalid_argument("ScenarioReductionConfig: CFLConfig cannot be null");
    }
    
    if (check_required) {
        // Check for mandatory members
        if (!config->has_member("k")) {
            throw std::invalid_argument("ScenarioReductionConfig: CFLConfig missing mandatory 'k' member");
        }
        
        // Validate k
        auto* k_config = dynamic_cast<const SimpleConfiguration<int>*>(config->get_member("k"));
        if (!k_config) {
            throw std::invalid_argument("ScenarioReductionConfig: k has wrong type in CFLConfig");
        }
        validate_k(k_config->get_value());
    }
    
    // Check optional members
    if (config->has_member("ell")) {
        auto* ell_config = dynamic_cast<const SimpleConfiguration<float>*>(config->get_member("ell"));
        if (!ell_config) {
            throw std::invalid_argument("ScenarioReductionConfig: ell has wrong type in CFLConfig");
        }
        validate_ell(ell_config->get_value());
    }
}

/*--------------------------------------------------------------------------*/

void ScenarioReductionConfig::validate_solver_config(const BlockSolverConfig* config, bool check_required)
{
    if (!config) {
        throw std::invalid_argument("ScenarioReductionConfig: SolverConfig cannot be null");
    }
    
    if (check_required) {
        // Check for mandatory members
        if (!config->has_member("algorithm")) {
            throw std::invalid_argument("ScenarioReductionConfig: SolverConfig missing mandatory 'algorithm' member");
        }
        
        // Validate algorithm
        auto* algorithm_config = dynamic_cast<const SimpleConfiguration<std::string>*>(config->get_member("algorithm"));
        if (!algorithm_config) {
            throw std::invalid_argument("ScenarioReductionConfig: algorithm has wrong type in SolverConfig");
        }
        validate_algorithm(algorithm_config->get_value());
    }
    
    // Validate ComputeConfig if present
    if (config->has_ComputeConfig(0)) {
        auto* compute_config = config->get_ComputeConfig(0);
        std::string computeName = config->get_ComputeName(0);
        
        // Check if the algorithm is MILP, which requires a specific solver
        auto* algorithm_config = dynamic_cast<const SimpleConfiguration<std::string>*>(config->get_member("algorithm"));
        if (algorithm_config && algorithm_config->get_value() == "MILP") {
            if (computeName != "HiGHSMILPSolver") {
                throw std::invalid_argument("ScenarioReductionConfig: Algorithm 'MILP' requires HiGHSMILPSolver, got '" + 
                                           computeName + "' instead");
            }
        } else if (computeName != "ScenarioReductionSolver" && computeName != "HiGHSMILPSolver") {
            throw std::invalid_argument("ScenarioReductionConfig: Invalid solver type '" + computeName + 
                                       "'. Must be 'ScenarioReductionSolver' or 'HiGHSMILPSolver'");
        }
    }
}

/*--------------------------------------------------------------------------*/

template<typename T>
SimpleConfiguration<T>* ScenarioReductionConfig::get_or_create_cfl_member(const std::string& name, const T& default_value)
{
    // Ensure we have a valid CFLConfig
    BlockConfig* cfl_config = ensure_cfl_config();
    
    // If member exists, check type
    if (cfl_config->has_member(name)) {
        auto* member = dynamic_cast<SimpleConfiguration<T>*>(cfl_config->get_member(name));
        if (member) {
            // Member has correct type, return it
            return member;
        } else {
            // Member has wrong type, replace it
            report_message("Replacing incompatible " + name + " member in CFLConfig");
            delete cfl_config->get_member(name);
            auto* new_member = new SimpleConfiguration<T>(default_value);
            cfl_config->add_member(name, new_member);
            return new_member;
        }
    } else {
        // Member doesn't exist, create it
        auto* new_member = new SimpleConfiguration<T>(default_value);
        cfl_config->add_member(name, new_member);
        return new_member;
    }
}

/*--------------------------------------------------------------------------*/

template<typename T>
SimpleConfiguration<T>* ScenarioReductionConfig::get_or_create_solver_member(const std::string& name, const T& default_value)
{
    // Ensure we have a valid SolverConfig
    BlockSolverConfig* solver_config = ensure_solver_config();
    
    // If member exists, check type
    if (solver_config->has_member(name)) {
        auto* member = dynamic_cast<SimpleConfiguration<T>*>(solver_config->get_member(name));
        if (member) {
            // Member has correct type, return it
            return member;
        } else {
            // Member has wrong type, replace it
            report_message("Replacing incompatible " + name + " member in SolverConfig");
            delete solver_config->get_member(name);
            auto* new_member = new SimpleConfiguration<T>(default_value);
            solver_config->add_member(name, new_member);
            return new_member;
        }
    } else {
        // Member doesn't exist, create it
        auto* new_member = new SimpleConfiguration<T>(default_value);
        solver_config->add_member(name, new_member);
        return new_member;
    }
}

// Explicit template instantiations for the types we'll use
template SimpleConfiguration<int>* ScenarioReductionConfig::get_or_create_cfl_member(const std::string&, const int&);
template SimpleConfiguration<float>* ScenarioReductionConfig::get_or_create_cfl_member(const std::string&, const float&);
template SimpleConfiguration<std::string>* ScenarioReductionConfig::get_or_create_solver_member(const std::string&, const std::string&);

/*--------------------------------------------------------------------------*/
/*------------------------- End ScenarioReductionConfig.cpp ---------------*/
/*--------------------------------------------------------------------------*/