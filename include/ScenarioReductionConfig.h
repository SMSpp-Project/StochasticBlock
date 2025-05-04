/*--------------------------------------------------------------------------*/
/*------------------- File ScenarioReductionConfig.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the ScenarioReductionConfig class, derived from Configuration,
 * which is intended to configure scenario reduction algorithms.
 *
 * \author Benoît Tran \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Benoît Tran
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __ScenarioReductionConfig
 #define __ScenarioReductionConfig
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Configuration.h"
#include "BlockConfig.h"
#include "BlockSolverConfig.h"
#include "ScenarioReductionSolver.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/
/// namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*-------------------- CLASS ScenarioReductionConfig -----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Configuration for scenario reduction algorithms
/** This is a ScenarioReductionConfig class that contains settings for the scenario reduction algorithms,
 * including parameters and algorithm selection. It follows the SMS++ configuration pattern by using
 * nested configuration objects rather than storing parameters directly in member variables.
 * 
 * The configuration can be loaded from/saved to a netCDF file with the following structure:
 * 
 * ScenarioReductionConfig  (Group)
 * ├── CFLConfig           (Group)
 * │   ├── k               (Attribute) = [int] Number of scenarios to select
 * │   └── ell             (Attribute) = [float] Power for Wasserstein distance (default: 2.0)
 * └── SolverConfig        (Group)
 *     └── algorithm       (Attribute) = [string] One of: "Dupacova" (default), "BestFit", "FirstFit", etc.
 * 
 * Explanation of parameters:
 * - k: Number of scenarios to select (must be positive, default: 10)
 * - ell: Power parameter in the ell-Wasserstein distance (must be positive, default: 2.0)
 * - algorithm: Method used for scenario reduction (default: "Dupacova")
 *   - "Baseline": Simple baseline method that selects scenarios based on probabilities
 *   - "Dupacova": Forward selection algorithm for discrete scenario reduction
 *   - "BestFit": Local search algorithm that selects best improvement at each step
 *   - "FirstFit": Local search algorithm that selects first satisfactory improvement
 *   - "MILP": Mixed Integer Linear Programming approach using HiGHSMILPSolver
 * 
 * Design Features:
 * - Nested Configuration Pattern: Uses two internal configuration objects:
 *   - f_cfl_config (BlockConfig): Stores parameters for CapacitatedFacilityLocationBlock (k, ell)
 *   - f_solver_config (BlockSolverConfig): Stores algorithm choice and solver settings
 * - Robust Parameter Management:
 *   - Helper methods (get_or_create_*) handle parameter access and creation
 *   - Validation methods ensure parameter values are valid
 *   - Default values used when parameters are missing or invalid
 * - Error Handling:
 *   - Centralized error reporting with consistent messages
 *   - Validation methods check parameter values
 *   - check_consistency() method verifies overall configuration integrity
 * - Consistency Management:
 *   - clone() method creates consistent deep copies
 *   - reset() method restores default configuration
 *   - auto-fix option in check_consistency() repairs invalid configurations
 * - Extraction Pattern:
 *   - extract_block_config() and extract_solver_config() methods generate specialized
 *     configuration objects for blocks and solvers based on stored settings
 *
 * Default Values:
 * - k = 10 (number of scenarios)
 * - ell = 2.0 (power parameter for Wasserstein distance)
 * - algorithm = "Dupacova" (scenario reduction method)
 */

class ScenarioReductionConfig : public Configuration
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------- CONSTRUCTING AND DESTRUCTING Configuration ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing ScenarioReductionConfig
 *  @{ */

 /// Default constructor (requires k parameter)
 /** Creates a ScenarioReductionConfig with the specified number of scenarios to select.
  *  Internally, creates a CFLConfig with the given k value and default ell=2.0,
  *  and a SolverConfig with default algorithm="Dupacova".
  *  @param k Number of scenarios to select (must be positive)
  */
 explicit ScenarioReductionConfig(int k);

/*--------------------------------------------------------------------------*/
 /// Construct from netCDF group
 /** Constructs a ScenarioReductionConfig from a netCDF group.
  * The group is expected to contain two mandatory subgroups:
  * - CFLConfig: Contains 'k' (mandatory) and 'ell' (optional) attributes
  * - SolverConfig: Contains 'algorithm' (mandatory) attribute or a full BlockSolverConfig
  * 
  * @param group The netCDF group containing scenario reduction settings
  * @throws std::invalid_argument If mandatory groups or attributes are missing
  */
 explicit ScenarioReductionConfig(const netCDF::NcGroup& group);

/*--------------------------------------------------------------------------*/
 /// Construct from parameters
 /** Constructs a ScenarioReductionConfig with the specified parameters.
  * Internally, creates a CFLConfig with the given k and ell values,
  * and a SolverConfig with the specified algorithm.
  * 
  * @param k Number of scenarios to select (must be positive)
  * @param ell Power parameter for Wasserstein distance (must be positive, default: 2.0)
  * @param algorithm Method to use (must be one of: "Baseline", "Dupacova", "BestFit", "FirstFit", "MILP", default: "Dupacova")
  */
 ScenarioReductionConfig(int k, float ell = 2.0f, const std::string& algorithm = "Dupacova");

/*--------------------------------------------------------------------------*/
 /// Clone method
 /** Creates a deep copy of this ScenarioReductionConfig instance.
  * @return A pointer to a new instance that is a copy of this one
  */
 [[nodiscard]] ScenarioReductionConfig* clone() const override;

/*--------------------------------------------------------------------------*/
 /// Destructor
 ~ScenarioReductionConfig() override;

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Configuration serialization and deserialization
 *  @{ */

 /// Deserialize from a netCDF group 
 /** Loads the configuration from a netCDF group structure.
  * The group is expected to contain two mandatory subgroups:
  * - CFLConfig: Contains 'k' (mandatory) and 'ell' (optional) attributes
  * - SolverConfig: Contains 'algorithm' (mandatory) attribute or a full BlockSolverConfig
  * 
  * @param group The netCDF group containing scenario reduction settings
  * @throws std::invalid_argument If mandatory groups or attributes are missing
  */
 void deserialize(const netCDF::NcGroup& group) override;

/*--------------------------------------------------------------------------*/
 /// Serialize to a netCDF group
 /** Saves the configuration to a netCDF group structure with two subgroups:
  * - CFLConfig: Contains 'k' and 'ell' attributes, plus any other members in f_cfl_config
  * - SolverConfig: Contains 'algorithm' attribute and a serialized f_solver_config
  * 
  * @param group The netCDF group to write settings to
  */
 void serialize(netCDF::NcGroup& group) const override;

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Accessors and mutators
 *  @{ */

 /// Set the number of scenarios to select in the CFLConfig
 /** Updates the k value in the internal CFLConfig. Creates a new CFLConfig if none exists.
  * @param value Number of scenarios (k) - must be positive
  * @throws std::invalid_argument If value is not positive
  */
 void set_k(int value);

 /// Get the number of scenarios to select from the CFLConfig
 /** @return Number of scenarios (k), or default value (10) if not found
  * @warning Prints a warning to std::cerr if using the default value
  */
 [[nodiscard]] int get_k() const;

/*--------------------------------------------------------------------------*/
 /// Set the power parameter for Wasserstein distance in the CFLConfig
 /** Updates the ell value in the internal CFLConfig. Creates a new CFLConfig if none exists.
  * @param value Power parameter (ell) - must be positive
  * @throws std::invalid_argument If value is not positive
  */
 void set_ell(float value);

 /// Get the power parameter for Wasserstein distance from the CFLConfig
 /** @return Power parameter (ell), or default value (2.0) if not found
  * @warning Prints a warning to std::cerr if using the default value
  */
 [[nodiscard]] float get_ell() const;

/*--------------------------------------------------------------------------*/
 /// Set the algorithm to use for scenario reduction in the SolverConfig
 /** Updates the algorithm value in the internal SolverConfig. Creates a new SolverConfig if none exists.
  * @param alg Algorithm name - must be one of: "Baseline", "Dupacova", "BestFit", "FirstFit", "MILP"
  * @throws std::invalid_argument If algorithm is not valid
  */
 void set_algorithm(const std::string& alg);

 /// Get the algorithm to use for scenario reduction from the SolverConfig
 /** @return Algorithm name, or default value ("Dupacova") if not found
  * @warning Prints a warning to std::cerr if using the default value
  */
 [[nodiscard]] const std::string& get_algorithm() const;

/*--------------------------------------------------------------------------*/
 /// Set the solver configuration
 /** Replaces the current SolverConfig with the provided one. If the new config
  * doesn't have an algorithm member, a default one will be added.
  * @param config BlockSolverConfig for the solver (ownership transferred)
  * @warning Prints a warning to std::cerr if adding a default algorithm
  */
 void set_solver_config(BlockSolverConfig* config);

 /// Get the solver configuration
 /** @return BlockSolverConfig for the solver (nullptr if not set)
  * @warning Prints a warning to std::cerr if returning nullptr
  */
 [[nodiscard]] BlockSolverConfig* get_solver_config() const;

 /// Reset configuration to default values
 /** Resets all configuration parameters to their default values.
  * This method:
  * - Deletes any existing configuration objects
  * - Creates new configuration objects with default values:
  *   - k = DEFAULT_K (10)
  *   - ell = DEFAULT_ELL (2.0)
  *   - algorithm = DEFAULT_ALGORITHM ("Dupacova")
  * - Reports a message to indicate the reset has been performed
  * 
  * Use this method to return to a known good state or to start fresh configuration.
  */
 void reset();
 
 /// Check consistency of the configuration
 /** Performs a thorough validation of the entire configuration structure.
  * Checks for:
  * - Missing or null configuration objects
  * - Missing required parameters
  * - Parameters with incorrect types
  * - Invalid parameter values
  * - Inconsistencies between parameters (e.g., MILP algorithm but wrong solver)
  * 
  * If auto_fix is true, attempts to repair any issues by:
  * - Creating missing configuration objects
  * - Adding missing parameters with default values
  * - Replacing parameters with wrong types
  * - Correcting invalid parameter values
  * - Adjusting inconsistent configurations
  * 
  * Reports all issues found, whether fixed or not.
  * 
  * @param auto_fix Whether to automatically fix inconsistencies (default: false)
  * @return true if the configuration is consistent (possibly after fixing), false otherwise
  */
 bool check_consistency(bool auto_fix = false);
 
/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Configuration extraction methods
 *  @{ */

 /// Extract a BlockConfig for CapacitatedFacilityLocationBlock
 /** Creates a BlockConfig for configuring a CFL block based on the scenario
  * reduction settings in this configuration. If f_cfl_config exists, it will be cloned
  * and extended with the required parameters. Otherwise, a new BlockConfig will be created.
  * 
  * The returned BlockConfig will have the following members:
  * - k: Number of scenarios to select (from f_cfl_config or get_k())
  * - ell: Power parameter for Wasserstein distance (from f_cfl_config or get_ell())
  * - n: Total number of scenarios (from the n_scenarios parameter)
  * 
  * @param n_scenarios Total number of scenarios in the original set
  * @return A BlockConfig for configuring a CFL block (caller owns the object)
  * @warning Prints warnings to std::cerr if using default values or replacing incompatible members
  */
 [[nodiscard]] BlockConfig* extract_block_config(size_t n_scenarios) const;

/*--------------------------------------------------------------------------*/
 /// Extract a BlockSolverConfig for scenario reduction
 /** Creates a BlockSolverConfig for configuring a solver based on the scenario
  * reduction settings in this configuration. If f_solver_config exists, it will be cloned
  * and extended with the required parameters. Otherwise, a new BlockSolverConfig will be created.
  * 
  * The returned BlockSolverConfig will:
  * - For algorithm="MILP": Configure a HiGHSMILPSolver with appropriate parameters
  * - For other algorithms: Configure a ScenarioReductionSolver with the specified algorithm
  * 
  * @return A BlockSolverConfig (caller owns the object)
  * @warning Prints warnings to std::cerr if using default values or making corrections
  */
 [[nodiscard]] BlockSolverConfig* extract_solver_config() const;

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*---------------------------- PROTECTED METHODS --------------------------*/

 /// Print the configuration to an output stream
 /** @param output The output stream to write to */
 void print(std::ostream& output) const override;

/*--------------------------------------------------------------------------*/
 /// Load from an input stream
 /** @param input The input stream to read from */
 void load(std::istream& input) override;

/*--------------------------------------------------------------------------*/
 /// Get list of valid algorithms for scenario reduction
 /** Returns a static vector containing all valid algorithm names from ScenarioReductionSolver::Algorithm.
  * Used for validating algorithm choices in various methods.
  * 
  * @return Vector of valid algorithm names: "Baseline", "Dupacova", "BestFit", "FirstFit", "MILP"
  */
 [[nodiscard]] static const std::vector<std::string>& get_valid_algorithms();
 
 /// Check if the given algorithm name is valid
 /** Validates that the algorithm name is one of the allowed values.
  * 
  * @param algorithm Algorithm name to validate
  * @return true if the algorithm is valid, false otherwise
  */
 [[nodiscard]] static bool is_valid_algorithm(const std::string& algorithm);
 
 /// Get a formatted string listing all valid algorithms
 /** @return A comma-separated list of valid algorithm names
  */
 [[nodiscard]] static std::string get_valid_algorithm_options();

/*--------------------------- PROTECTED FIELDS  ---------------------------*/

 BlockConfig* f_cfl_config;        ///< CFLConfig: contains k and ell parameters
 BlockSolverConfig* f_solver_config; ///< SolverConfig: contains algorithm and solver settings

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 private:
 
 /// Default values for configuration parameters
 static constexpr int DEFAULT_K = 10;              ///< Default number of scenarios
 static constexpr float DEFAULT_ELL = 2.0f;        ///< Default power parameter
 static const inline std::string DEFAULT_ALGORITHM = "Dupacova"; ///< Default algorithm

 /// Initialize the configuration with default values
 /** Common initialization method called by all constructors to ensure consistent
  * initialization of the configuration. Creates empty configuration objects
  * with default values if they don't exist yet.
  * 
  * @param k Initial value for k parameter (must be positive)
  * @param ell Initial value for ell parameter (must be positive)
  * @param algorithm Initial value for algorithm (must be valid)
  * @throws std::invalid_argument If parameters are invalid
  */
 void initialize(int k, float ell, const std::string& algorithm);

 /// Helper method to ensure CFLConfig exists and is valid
 /** Ensures that f_cfl_config is not null, creating it if needed.
  *
  * @return A pointer to the CFLConfig object (never null)
  */
 BlockConfig* ensure_cfl_config();

 /// Helper method to ensure SolverConfig exists and is valid
 /** Ensures that f_solver_config is not null, creating it if needed.
  *
  * @return A pointer to the SolverConfig object (never null)
  */
 BlockSolverConfig* ensure_solver_config();

 /// Helper method to get or create a member in the CFLConfig
 /** Gets an existing member from the CFLConfig or creates a new one if it doesn't exist.
  * If the member exists but has the wrong type, it will be replaced.
  *
  * @tparam T Type of the configuration value
  * @param name Name of the member
  * @param default_value Default value to use if the member needs to be created
  * @return A pointer to the SimpleConfiguration object (never null)
  */
 template<typename T>
 SimpleConfiguration<T>* get_or_create_cfl_member(const std::string& name, const T& default_value);

 /// Helper method to get or create a member in the SolverConfig
 /** Gets an existing member from the SolverConfig or creates a new one if it doesn't exist.
  * If the member exists but has the wrong type, it will be replaced.
  *
  * @tparam T Type of the configuration value
  * @param name Name of the member
  * @param default_value Default value to use if the member needs to be created
  * @return A pointer to the SimpleConfiguration object (never null)
  */
 template<typename T>
 SimpleConfiguration<T>* get_or_create_solver_member(const std::string& name, const T& default_value);
 
 /// Method for consistent error reporting
 /** Used for reporting warnings and errors with consistent formatting.
  * Error messages are prefixed with the class name for easier identification.
  * 
  * @param message Message to report
  * @param severity Severity level: 0 = info, 1 = warning, 2 = error
  */
 static void report_message(const std::string& message, int severity = 1);
 
 /// Helper method for reporting when using default values
 /** Reports a consistent warning message when using a default value.
  * 
  * @param param_name Name of the parameter
  * @param default_value Default value being used
  */
 static void report_using_default(const std::string& param_name, const std::string& default_value);
 
 /// Validate the k parameter (number of scenarios)
 /** Validates that k is a positive integer. This method performs a simple
  * validation check to ensure that k > 0. Used by accessors, mutators,
  * and consistency checking to maintain valid configuration.
  * 
  * @param k Value to validate
  * @throws std::invalid_argument If k is not positive (k <= 0)
  */
 static void validate_k(int k);
 
 /// Validate the ell parameter (power for Wasserstein distance)
 /** Validates that ell is a positive float. This method performs a simple
  * validation check to ensure that ell > 0.0. Used by accessors, mutators,
  * and consistency checking to maintain valid configuration.
  * 
  * @param ell Value to validate
  * @throws std::invalid_argument If ell is not positive (ell <= 0.0)
  */
 static void validate_ell(float ell);
 
 /// Validate the algorithm parameter
 /** Validates that the algorithm is one of the supported values listed in 
  * get_valid_algorithms(). This method uses is_valid_algorithm() internally
  * to check whether the given algorithm is supported.
  * 
  * @param algorithm Value to validate
  * @throws std::invalid_argument If algorithm is not in the list of valid algorithms
  * @see get_valid_algorithms(), get_valid_algorithm_options(), is_valid_algorithm()
  */
 static void validate_algorithm(const std::string& algorithm);
 
 /// Validate a CFLConfig
 /** Performs comprehensive validation of a CFLConfig object.
  * Validates that:
  * - The config pointer is not null
  * - If check_required is true, checks for mandatory members (k)
  * - Validates that member types are correct (e.g., k is int, ell is float)
  * - Validates that member values are valid (e.g., k > 0, ell > 0)
  * 
  * @param config BlockConfig to validate
  * @param check_required Whether to check for required members (default: true)
  * @throws std::invalid_argument If the config is invalid
  */
 static void validate_cfl_config(const BlockConfig* config, bool check_required = true);

 /// Validate a SolverConfig
 /** Performs comprehensive validation of a SolverConfig object.
  * Validates that:
  * - The config pointer is not null
  * - If check_required is true, checks for mandatory members (algorithm)
  * - Validates that member types are correct
  * - Validates that member values are valid (algorithm is supported)
  * - Checks consistency between algorithm and ComputeConfig if present
  *   (e.g., MILP algorithm requires HiGHSMILPSolver)
  * 
  * @param config BlockSolverConfig to validate
  * @param check_required Whether to check for required members (default: true) 
  * @throws std::invalid_argument If the config is invalid
  */
 static void validate_solver_config(const BlockSolverConfig* config, bool check_required = true);
 
 /// Helper method to deserialize the CFLConfig part
 /** Extracts and parses the CFLConfig subgroup from a netCDF group.
  * This helper method is used by the main deserialize() method to modularize
  * the deserialization process. It follows these steps:
  * 1. Locates the CFLConfig subgroup in the provided netCDF group
  * 2. Extracts and validates the 'k' parameter (mandatory)
  * 3. Extracts and validates the 'ell' parameter (optional, defaults to DEFAULT_ELL)
  * 4. Creates a new BlockConfig with these parameters
  * 
  * @param group The netCDF group containing the ScenarioReductionConfig
  * @return A pointer to the newly created BlockConfig (caller owns the object)
  * @throws std::invalid_argument If the CFLConfig subgroup is missing or required data is invalid
  */
 BlockConfig* deserialize_cfl_config(const netCDF::NcGroup& group) const;
 
 /// Helper method to deserialize the SolverConfig part
 /** Extracts and parses the SolverConfig subgroup from a netCDF group.
  * This helper method is used by the main deserialize() method to modularize
  * the deserialization process. It handles two formats:
  * 1. A full BlockSolverConfig object (with 'type' attribute)
  * 2. A simpler configuration with just an 'algorithm' attribute
  * 
  * For the first case, it uses the factory pattern to create the appropriate
  * configuration object. For the second case, it creates a new BlockSolverConfig
  * with just the algorithm parameter.
  * 
  * @param group The netCDF group containing the ScenarioReductionConfig
  * @return A pointer to the newly created BlockSolverConfig (caller owns the object)
  * @throws std::invalid_argument If the SolverConfig subgroup is missing or required data is invalid
  */
 BlockSolverConfig* deserialize_solver_config(const netCDF::NcGroup& group) const;

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

};  // end( class( ScenarioReductionConfig ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- FACTORY MANAGEMENT ----------------------------*/
/*--------------------------------------------------------------------------*/

// Register ScenarioReductionConfig to the Configuration factory
// (This will be defined in the .cpp file)
// SMSpp_insert_in_factory_cpp_0( ScenarioReductionConfig );

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* ScenarioReductionConfig.h included */

/*--------------------------------------------------------------------------*/
/*----------------- End File ScenarioReductionConfig.h ---------------------*/
/*--------------------------------------------------------------------------*/