/*--------------------------------------------------------------------------*/
/*--------------------- File ScenarioGenerator.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the definition of the *abstract* classes ScenarioGenerator
 * and MultiStageScenarioGenerator, that represent the interface between a
 * user that needs to generate instantiations of random processes (scenarios),
 * typically for the purpose of forming / solving stochastic optimization
 * problems, and the actual data structures holding the data that contains /
 * generates these instantiations.
 *
 * There can be several different ways in which such data is organized,
 * starting from the fact that there is only one "stage" (all the random
 * variables are revealed at the same time) or multiple ones (there is a
 * different random variable for each of a discrete set of points in time,
 * called stages, and random variables are progressively revealed as time
 * passes through them). Once the point in time is fixed, though, the data
 * looks the same in both cases; this is why the base class ScenarioGenerator
 * defines the interface for retrieving the data concerning a single stage,
 * and MultiStageScenarioGenerator derives from it by adding the concepts
 * related to handling the different points in time (stages).
 *
 * Apart from that, it is expected that different concrete derived classes
 * will be implemented to consider relevant special cases, such as:
 *
 * - DiscreteScenarioSet deriving from ScenarioGenerator for the case where
 *   the (single-stage) scenarios are just the result of a simple unique
 *   random variable with a discrete distribution (one element for each
 *   scenario, just provided in input in their final form).
 *
 * - MultiStageDiscreteScenarioSet deriving from MultiStageScenarioGenerator
 *   for the case where the (multistage) scenarios are just the result of a
 *   simple unique random variable with a discrete distribution (one element
 *   for each scenario, just provided in input in their final form).
 *
 * - IndependentMultiStageScenarioGenerator where each stage is a random
 *   variable independent from the past history; a possible implementation is
 *   just as a vector of ScenarioGenerator *, one for each of the stages.
 *
 * Many other possibilities exist, e.g., single-stage scenarios constructed
 * from the combination of an arbitrary number of distributions (discrete,
 * continuous of various types, ...), multi-stage scenarios in the Markovian
 * setting where each stage is (an arbitrary combination of random variables)
 * depending only from the results of the previous stage, general multi-stage
 * scenarios where the random variables depend on the entire previous history,
 * and so on. In many cases the user is not concerned with these details, and
 * therefore they can just program against the [MultiStage]ScenarioGenerator
 * interface.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Benoit Tran \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni and Benoit Tran
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __ScenarioGenerator
 #define __ScenarioGenerator
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SMSTypedefs.h"
#include "Configuration.h" // For Configuration class

// Required for Scenario type but not available in SMSTypedefs.h yet
// TODO: Consider adding to SMSTypedefs.h in a future update
#include <memory>    // for std::unique_ptr (View handles)
#include <span>      // C++20: For non-owning views of contiguous data
#include <stdexcept> // for std::logic_error / std::invalid_argument
#include <vector>    // for std::vector< ScenarioIndex > overloads

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
// Forward declaration of Block class
class Block;

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASSES -----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Classes in ScenarioGenerator
 *  @{ */

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS ScenarioGenerator ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// abstract base class for a handler of scenarios
/** ScenarioGenerator define an abstract interface that separates "users of
 * scenarios", i.e., components which need to access specific realizations of
 * (sets of) random variables, from the implementation of the random variables
 * themselves. It is designed to work in tandem with its derived (still
 * abstract) class MultiStageScenarioGenerator, hence these comments actually
 * address in part both. This is because both share the interface that allows
 * to access the data (realization) of one specific point in time (stage), of
 * which there is only one in the case of ScenarioGenerator.
 *
 * In the general framework of a multistage problem, a number T > 0 of time
 * instant (stages) are defined. Then, a scenario ( x_0 , ... , x_{T-1} ) is a
 * realization of a finite horizon stochastic process ( X_0 , ... , X_{T-1} ),
 * where for each 0 <= t <= T - 1, the random variable X_t has values in some
 * Euclidean space R^(d_t). The base class makes no assumption on the
 * "internal structure" of x_t, which can in general be the composition of
 * different random variables (say, a pair ( c , a ) where c is the vector
 * c[ j ] denoting the cost of a given commodity at a predetermined set of
 * locations j \in J, and a is the vector a[ j ] denoting the availability of
 * the same commodity at the same places; or the matrix of pairs where each
 * row correspond to a different commodity, ...). As far as the base classes
 * are concerned, a scenario x_t for a given stage t is just a vector of
 * appropriate length d_t, whose internal structure is supposed to be fixed
 * and known to the user. Also, the base classes make no assumption on how
 * the variable X_t depends on ( x_0 , ... , x_{t-1} ) (the history at stage
 * t).
 *
 * When T = 1, the time dimension is implicit, i.e., the stochastic process is
 * its only element, the random variable X_0, and a scenario is its unique
 * value x_0. This is why there are two separate classes: ScenarioGenerator
 * and MultiStageScenarioGenerator. In the former the time dimension is
 * fixed, and therefore there is no need of methods to manage it, differently
 * from the latter. However, once the stage t is known in the multistage
 * setting then a scenario x_t for that stage "looks a lot like" in the
 * single-stage case, which is why a large part of the interface is shared.
 *
 * Other notable elements of ScenarioGenerator:
 *
 * - it has a factory (like, Block, Solver, Configuration, ...) managed by
 *   means of the standard SMSpp_insert_in_factory_h and
 *   SMSpp_insert_in_factory_cpp_[0_t]() macros;
 *
 * - it therefore has a static method new_ScenarioGenerator( std::string )
 *   that constructs a :ScenarioGenerator of the specific type dictated by
 *   the string argument (using the factory)
 *
 * - it has an in-built hierarchy of deserialization methods, starting with
 *   deserialize( std::string ) that deserializes a :ScenarioGenerator from
 *   the netCDF file with the given filename, passing through
 *   new_ScenarioGenerator( netCDF::NcGroup ) that reads some top-level
 *   information from the group (either the class type or a filename),
 *   finally down to the pure virtual deserialize( netCDF::NcGroup ) that
 *   actually implements the deserialization for the specific
 *   :ScenarioGenerator
 *
 * - it allows to set the seed of the pseudo-random number generator(s)
 *   possibly involved in the (sampling) operations
 *
 * - distinguishes different types of "scenario pools" supporting different
 *   use cases
 *
 * - ... ?
 */

class ScenarioGenerator
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *  @{ */

/// type for reading the data characterizing a scenario.
/** Type for reading the data characterizing a scenario. The base class does
 * not make any assumption on any "structure" that a scenario may have; a
 * scenario is just a vector of numbers (double). It is assumed that the
 * representation of the "current scenario" that can be read as a contiguous
 * array of double will always be available inside any implementation of the
 * abstract base class, hence a Scenario is just a std::span<> provided
 * read-only access to that data structure. */

 using Scenario = std::span< const double >;

/*--------------------------------------------------------------------------*/
/// type for the index characterizing a specific scenario
/** Integer type that can be used to identify one of the possible finite
 * number of realizations that can be produced. Basically implies the
 * maximum number of them. */

 using ScenarioIndex = unsigned int;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// constexpr for "infinitely many scenarios"
 static constexpr ScenarioIndex INFScenario = Inf< ScenarioIndex >();

/*--------------------------------------------------------------------------*/
/// type for the scenario size (the length of the double vector)
/** Integer type that can be used to specify the length of the double vector
 * containing the data of the scenario. */

 using ScenarioSize = size_t;

/*--------------------------------------------------------------------------*/
/// Note about Block usage
/** Derived classes of ScenarioGenerator *might* need the data of a problem
 * contained into a Block. These derived class will include (derived) classes
 * of Block. The Block class is forward declared at the namespace level. */

/** @} ---------------------------------------------------------------------*/
/*------------ CONSTRUCTING AND DESTRUCTING ScenarioGenerator --------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing ScenarioGenerator
 *  @{ */

 /// constructor: initialize the few data structures of the base class
 ScenarioGenerator( void ) { }

/*--------------------------------------------------------------------------*/
 /// construct a :ScenarioGenerator of a specific type using the factory
 /** Use the ScenarioGenerator factory to construct a :ScenarioGenerator
  * object of a type specified by classname (a std::string with the name of
  * the class inside). Note that the method is static because the factory is
  * static, hence it is to be called as
  *
  *   ScenarioGenerator *mySG = ScenarioGenerator::new_ScenarioGenerator(
  *                                                             some_class );
  *
  * i.e., without any reference to any specific ScenarioGenerator.
  *
  * For this to work, each :ScenarioGenerator has to:
  *
  * - add the line
  *
  *     SMSpp_insert_in_factory_h;
  *
  *   to its definition (typically, in the private part in its .h file);
  *
  * - add the line
  *
  *     SMSpp_insert_in_factory_cpp_0( name_of_the_class );
  *
  *   to exactly *one* .cpp file, typically that :ScenarioGenerator .cpp file.
  *   If the name of the class contains any parentheses, then one must enclose
  *   the name of the class in parentheses and instead add the line
  *
  *     SMSpp_insert_in_factory_cpp_0( ( name_of_the_class ) );
  *
  * Any whitespaces that the given \p classname may contain is ignored.
  *
  * @param classname The name of the specific :ScenarioGenerator class that
  * must be constructed. */

 [[nodiscard]] static ScenarioGenerator * new_ScenarioGenerator(
					   const std::string & classname ) {
  const std::string classname_( SMSpp_classname_normalise(
					        std::string( classname ) ) );
  const auto it = ScenarioGenerator::f_factory().find( classname_ );
  if( it == ScenarioGenerator::f_factory().end() )
   throw( std::invalid_argument( classname +
			     " not present in ScenarioGenerator factory" ) );
  return( ( it->second )() );
  }

/*--------------------------------------------------------------------------*/
 /// de-serialize a :ScenarioGenerator out of a file
 /** Top-level de-serialization method: takes the \p filename of a netCDF file
  * and returns the complete :ScenarioGenerator object whose description is
  * found in the file. Uses the new_ScenarioGenerator( netCDF::NcGroup )
  * second-level method to allow indirect file access, see the comment there
  * for the format of the file.
  *
  * If anything goes wrong with deserialization, nullptr is returned.
  *
  * Note that the method is static, hence it is to be called as
  *
  *      auto mySG = ScenarioGenerator::deserialize( some_file );
  *
  * i.e., without any reference to any specific ScenarioGenerator. */

 [[nodiscard]] static ScenarioGenerator * deserialize(
					     const std::string & filename) {
  try {
   // Relies on netCDF API to handle is_open() check
   /* might need to change filename into a cstring, if it works,
    * let antonio knows */
   netCDF::NcFile dataFile( filename , netCDF::NcFile::read );

   return( ScenarioGenerator::new_ScenarioGenerator( dataFile ) );
   }
  catch( netCDF::exceptions::NcException & e ) {
   std::cerr << "Error opening or processing netCDF file: " << e.what()
	     << std::endl;
   }
  catch( std::exception & e ) {
   std::cerr << "error " << e.what() << " in deserialize" << std::endl;
   }
  catch( ... ) {
   std::cerr << "unknown error in deserialize" << std::endl;
   }
  return( nullptr );
  }

/*--------------------------------------------------------------------------*/
 /// de-serialize a :ScenarioGenerator out of an netCDF::NcGroup, returns it
 /** Second-level de-serialization method: takes a netCDF::NcGroup supposedly
  * containing  (all the information describing) a :ScenarioGenerator (either
  * "directly" or "indirectly") and returns a pointer to a newly minted
  * :ScenarioGenerator object corresponding to what is found in the file.
  *
  * The method works with two different kinds of netCDF::NcGroup:
  *
  * - A "direct" group that contains at least the string attribute "type";
  *   this is used in the factory to construct an "empty" :ScenarioGenerator
  *   of that type [see new_ScenarioGenerator( std::string & )], and then the
  *   method deserialize( netCDF::NcGroup ) of the newly minted
  *   :ScenarioGenerator is invoked (with argument \p group) to finish the
  *   work.
  *
  * - An "indirect" group that just needs to contain the single string
  *   attribute "filename"; in this case, the attribute is used as argument
  *   for a call to deserialize( const std::string & ) that will extract the
  *   :ScenarioGenerator by the corresponding netCDF file.
  *
  * In the case \p group contains both "type" and "filename", the first takes
  * the precedence (direct groups have precedence over indirect ones).
  *
  * Note that this method is static (see the previous versions for comments
  * about it) and returns a pointer to ScenarioGenerator, hence it has to
  * have a different name from deserialize( netCDF::NcGroup ) (since the
  * signature is the same but for the return type).
  *
  * If anything goes wrong with the process, nullptr is returned. */

 [[nodiscard]] static ScenarioGenerator * new_ScenarioGenerator(
					     const netCDF::NcGroup & group )
 {
  try {
   if( group.isNull() )
    return( nullptr );

   std::string tmp;
   auto gtype = group.getAtt( "type" );
   if( gtype.isNull() ) {
    auto gfile = group.getAtt( "filename" );
    if( gfile.isNull() )
     return( nullptr );

    gfile.getValues( tmp );

    return( deserialize( tmp ) );
    }

   gtype.getValues( tmp );
   auto result = new_ScenarioGenerator( tmp );
   result->deserialize( group );
   return( result );
   }
  catch( netCDF::exceptions::NcException & e ) {
   std::cerr << "netCDF error " << e.what() << " in deserialize" << std::endl;
   }
  catch( std::exception & e ) {
   std::cerr << "error " << e.what() << " in deserialize" << std::endl;
   }
  catch( ... ) {
   std::cerr << "unknown error in deserialize" << std::endl;
   }

  return( nullptr );
  }

/*--------------------------------------------------------------------------*/
 /// de-serialize the current :ScenarioGenerator out of netCDF::NcGroup
 /** Third and final level de-serialization method: takes a netCDF::NcGroup
  * supposedly containing all the information required to de-serialize the
  * ScenarioGenerator, and initialize the current ScenarioGenerator out of it.
  *
  * A group containing a :ScenarioGenerator must have a mandatory string
  * attribute "type" that contains the classname() of the
  * :ScenarioGenerator, which is actually useful at the higher level of the
  * deserialize() hierarchy where the :ScenarioGenerator has yet to be
  * constructed, rather than at this point where it clearly already has,
  * plus whatever other information is required by the specific
  * :ScenarioGenerator.
  *
  *      THIS IS THE METHOD TO BE IMPLEMENTED BY DERIVED CLASSES
  *
  * and in fact it is pure virtual. */

 virtual void deserialize( const netCDF::NcGroup & group ) = 0;

/*--------------------------------------------------------------------------*/
 /// serialize the current ScenarioGenerator to netCDF::NcGroup
 /** Method to serialize the ScenarioGenerator to a netCDF::NcGroup.
  * This is the counterpart to deserialize() and should save all the
  * information required to reconstruct the ScenarioGenerator state.
  *
  * The base class implementation only writes the mandatory "type"
  * attribute (the classname(), needed by new_ScenarioGenerator() to
  * reconstruct the :ScenarioGenerator out of the group). Derived classes
  * should override this method to save their specific data, calling the
  * base class method first.
  *
  * @param group The netCDF group to serialize to */

 virtual void serialize( netCDF::NcGroup & group ) const {
  group.putAtt( "type" , classname() );
  }

/*--------------------------------------------------------------------------*/
 /// destructor

 virtual ~ScenarioGenerator() = default;

/** @} ---------------------------------------------------------------------*/
/*----- METHODS FOR READING THE STATIC DATA OF THE ScenarioGenerator -------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the static data of the ScenarioGenerator
 *  @{ */

 /// returns the maximum number of scenarios that can be generated
 /** This method returns the support size of the random variable of the
  * ScenarioGenerator, i.e., the maximum number of different scenarios that it
  * can ever output. This is done so that the user avoids calling the
  * init_*_pool() methods [see] with requiring too high a number of scenarios.
  * In the case of the random variable being continuous, the support size is
  * (theoretically) infinite, so INFScenario is reported. This is done by the
  * base class implementation, so that the method only needs to be redefined
  * by derived classes implementing variables with finite support. */

 [[nodiscard]] virtual ScenarioIndex get_support_size( void ) {
  return( INFScenario );
  }

/*--------------------------------------------------------------------------*/
 /// returns the size of any scenario (length of the double vector)
 /** This method returns the size of the double vector containing the data
  * that specify the instantiation x_0 of the random variable X_0. This will
  * be the size of the std::span< const double > (Scenario) returned by
  * get_current_scenario(). */

 [[nodiscard]] virtual ScenarioSize get_scenario_size( void ) const = 0;

/*--------------------------------------------------------------------------*/
 /// getting the classname of this ScenarioGenerator
 /** Given a :ScenarioGenerator (i.e., any class derived by it), this method
  * returns a string with its class name; unlike std::type_info.name(), there
  * *are* guarantees, i.e., the name will always be the same.
  *
  * The method works by dispatching the private virtual method private_name().
  * The latter is automatically implemented by the
  * SMSpp_insert_in_factory_cpp_* macros [see SMSTypedefs.h], hence this
  * comes at no cost since these have to be called somewhere to ensure that
  * any :ScenarioGenerator will be added to the factory. Actually, since
  * ScenarioGenerator::private_name() is pure virtual, this ensures that it
  * is not possible to forget to call SMSpp_insert_in_factory_cpp_* for any
  * :ScenarioGenerator because otherwise it is a pure virtual class (unless
  * the programmer purposely defines private_name() without calling the
  * macro, which seems rather pointless). */

 [[nodiscard]] const std::string & classname( void ) const {
  return( private_name() );
  }

/*--------------------------------------------------------------------------*/
 /// setting a partner Block
 /** Although ScenarioGenerator operations should be completely independent
  * on the specific type of :Block in which the scenario gets realised,
  * there may be cases in which some dependency may be necessary / useful.
  * For instance, in the context of scenario reduction (cf.
  * init_representative_pool() and the comments therein), one may want to
  * construct a problem-dependent metric to guide the scenario reduction
  * process. In such case, the definition of the tailored problem-dependent
  * metric depends on some of the problem's data that needs to be extracted
  * from the :Block. This means that 1) a specialised :ScenarioGenerator
  * must be written to handle this, and 2) it has to be provided with a
  * pointer to the original Block.
  *
  * ScenarioGenerator provides this virtual method for this purpose.
  * Whomever is building a ScenarioGenerator should pass it the pointer to
  * the interested Block via this method, that by default does nothing.
  * Derived classes that need information from the Block will have to
  * overwrite the method and check dynamically that the given Block pointer
  * is indeed if the expected type.
  *
  * A main user of ScenarioGenerator is StochasticBlock, so it will be
  * it providing ScenarioGenerator with the Block pointer. In this case,
  * the passed Block will be the StochasticBlock itself rather than its
  * inner Block. This allows access to both:
  * - the inner Block, via get_inner_block(), after the pointer is
  *   dynamic_cast< StochasticBlock * >-ed;
  * - the stochastic structure / scenarios managed by the StochasticBlock.
  *
  * Thus, the overriding method should:
  * 1. Check if the Block is a StochasticBlock (via dynamic_cast)
  * 2. Access the inner Block via get_inner_block()
  * 3. Verify the inner Block is of the expected type
  * 4. Verify that scenarios at ScenarioGenerator's disposal are coherent
  *    with the expected stochasticity of the Block's problem
  *
  * For instance, a StochasticBlock wrapping a
  * CapacitatedFacilityLocationBlock might have stochasticity in either
  * demands or capacities. The corresponding derived version of set_Block
  * should then take care of checking that the scenarios generated by the
  * derived class of ScenarioGenerator are coherent with demands or
  * capacities. */

  virtual void set_Block( Block * ) {};

/*--------------------------------------------------------------------------*/
 /// setting configuration for the ScenarioGenerator
 /** This method allows setting configuration parameters for the
  * ScenarioGenerator. The base class implementation does nothing, but
  * derived classes can override this to accept Configuration objects.
  *
  * The Configuration object passed should contain the appropriate
  * parameters for the specific :ScenarioGenerator implementation.
  * Implementations should validate the Configuration and may throw
  * exceptions if invalid parameters are provided.
  *
  * @param config The Configuration object containing parameters */

  virtual void set_config( Configuration* config ) {}

/** @} ---------------------------------------------------------------------*/
/*-------------------- METHODS FOR READING THE SCENARIOS -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Generating the scenario pool and reading the scenarios
 *
 * ScenarioGenerator works based on the concept of "scenario pool". The user
 * has to request the (logical) creation of a finite set of scenarios that
 * they can then read. Scenario pools can be constructed differently with
 * different aims in mind (and a different cost). Once a scenario pool is
 * initialized, it can be traversed by going through the scenarios in the
 * pool one by one; at each stage a current scenario is defined that can be
 * read. Once the current scenario is moved forward, the previous current
 * scenario is in principle no longer available.
 *
 * The user is supposedly incapable of differentiating the scenarios of the
 * pool from one another, which means that there is no assumption on the
 * order in which the scenarios of the pool are visited. Also, there should
 * be no need to re-visit a scenario once the user has visited it and moved
 * on, so "the clock of the pool only moves forward" (if the user needs to
 * look at more than one scenario at a time, they will have to save them in
 * their own data structure.
 *  @{ */

 /// set the seed of the internal random generator
 /** Creating pools can well involve random sampling operations, and therefore
  * the use of pseudo-random numbers. To guarantee reproducibility, the seed
  * for the pseudo-random generator can be set with this method. It is
  * conceivable that some "complicated" ScenarioGenerator may require more
  * than one generator, if this is the case, the method may have to be
  * changed (but it is not impossible to use the same seed for all of them,
  * or to use the seed within a random number generator that is used to
  * generate the seeds, so just one seed may well be enough). */

 virtual void set_seed( unsigned long seed = 0 ) = 0;

/*--------------------------------------------------------------------------*/
 /// generate a random pool by sampling from the current universe
 /** This method (logically) constructs the scenario pool by sampling \p
  * size scenarios from the *current universe* of the generator. The
  * universe is the (possibly restricted) set of scenarios available to
  * the random sampling:
  *
  *  - by default (no prior call to init_representative_pool()) the
  *    universe is the full set of scenarios the generator can produce;
  *
  *  - a previous call to init_representative_pool(K) restricts the
  *    universe to K representatives (the "filter"). The restriction is
  *    *persistent* across subsequent init_random_pool() calls;
  *
  *  - calling init_representative_pool(INFScenario) resets the universe
  *    back to the full set.
  *
  * The scenarios in the resulting pool are typically used during a
  * simulation phase to evaluate the quality of decisions taken by some
  * optimization model that has only "seen" a (much) smaller, more
  * carefully chosen set of representative scenarios.
  *
  * \p size is the number of scenarios in the produced pool. The
  * default value INFScenario means "use the whole current universe"
  * (i.e., size := |universe|). Any finite value of \p size must be
  * less than or equal to the size of the current universe; otherwise
  * an exception is thrown (it is a caller bug to ask for more
  * scenarios than the universe can provide — use INFScenario to opt
  * into the "give me everything" case explicitly).
  *
  * Once the method returns, the first scenario in the pool is available
  * as the current scenario and can be read [see get_current_scenario()
  * and get_current_scenario_probability()] right away.
  *
  * In the :MultiStageScenarioGenerator setting the stage is the one of
  * the View this is called on [see MultiStageScenarioGenerator::View],
  * consistently with the way every other "current scenario" method is
  * interpreted (get_scenario_size(), get_current_scenario(),
  * next_scenario(), ...). To (re-)init the pool of more than one stage,
  * the caller walks the tree with View::descend() and calls
  * init_random_pool() once per stage. Being a pool re-definition, this is
  * one of the two operations that *write* the process, and it invalidates
  * the Views that were relying on the pool it re-defines [see
  * MultiStageScenarioGenerator::View]. */

  virtual void init_random_pool( ScenarioIndex size = INFScenario ) = 0;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// restrict the universe to a "representative" subset of given size
 /** This method (logically) restricts the *universe* of available
  * scenarios to a subset of \p size representatives, working hard to
  * ensure that the chosen subset is "as much representative as
  * possible" of the whole (possibly, infinite) set of possible
  * realizations. This may be nontrivial: for the discrete case it
  * typically involves solving an Optimal Transport Problem, hence a
  * potentially non-negligible computational cost. The restriction is
  * persistent and survives across later calls to init_random_pool(),
  * which will then sample from these representatives only. To remove
  * the restriction and return to the full universe, call this method
  * again with \p size == INFScenario (the default).
  *
  * As a side effect, this method also initialises the current pool to
  * the universe just selected, walked in its canonical (i.e., natural)
  * order — there is no shuffling. Hence after a call to
  * init_representative_pool(K) the user can directly read the K
  * representatives with reset_pool() + next_scenario() in their
  * canonical order, without having to call init_random_pool()
  * separately. This is what the TwoStageStochasticBlock-style users
  * want: "give me all my scenarios in order".
  *
  * The decisions taken by a model built on the restricted universe can
  * subsequently be evaluated with init_random_pool(K') (with K' larger
  * than \p size and limited to the restricted universe) to draw a
  * random sample for simulation.
  *
  * \p size is the number of representatives. The default value
  * INFScenario means "the universe is the full set of scenarios"
  * (i.e., no restriction): use it to *reset* a previous restriction.
  * Other special values: a value larger than the natural support size
  * is treated as INFScenario.
  *
  * Once the method returns, the first scenario in the pool is
  * available as the current scenario and can be read [see
  * get_current_scenario() and get_current_scenario_probability()]
  * right away.
  *
  * In the :MultiStageScenarioGenerator setting,
  * init_representative_pool() acts on the pool of the View this is
  * called on; see init_random_pool() above for the rationale. To
  * (re-)init the representative universe of more than one stage, walk
  * the tree with View::descend() and call init_representative_pool() at
  * each. In a genuine scenario *tree* this is the very operation that
  * selects which children a node has, hence it has to be called at a node
  * before its children can be read, and it invalidates the Views that were
  * relying on the pool it re-defines [see MultiStageScenarioGenerator::
  * View]. A generator that materialises the whole tree upfront can instead
  * reduce it once and for all as a whole, in which case the call at a node
  * has nothing left to do; see the specific derived class. */

  virtual void init_representative_pool( ScenarioIndex size = INFScenario ) = 0;

/*--------------------------------------------------------------------------*/
 /// read the data of the current scenario
 /** Gets a std::span< const double > (Scenario) of size get_scenario_size()
  * [see] that contains the data that specify the instantiation x_0 of the
  * random variable X_0. This points to an internal data structure of the
  * :ScenarioGenerator, which is guaranteed to remain available up until the
  * next call of either next_scenario() or init_*_pool(), after which the data
  * structure may no longer be available and the span should no longer be
  * used. */

 [[nodiscard]] virtual Scenario get_current_scenario( void ) const = 0;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// return the probability associated to the current scenario

 [[nodiscard]] virtual double get_current_scenario_probability( void )
  const = 0;

/*--------------------------------------------------------------------------*/
 /// move the current scenario to the next scenario in the pool
 /** Tries to move forward the current scenario to the next scenario in the
  * pool, returns true if this succeeds, which should always happen if the
  * number of scenarios in the pool seen so far is lower than the size
  * specified in init_*_pool(). Upon receiving a true response, the user can
  * call again get_current_scenario[_probability]() to retrieve the data of
  * the new scenario, knowing that any Scenario corresponding to the previous
  * scenarios seen has potentially been invalidated and should no longer be
  * used. */

 virtual bool next_scenario( void ) = 0;

/*--------------------------------------------------------------------------*/
 /// reset the pool iteration to the beginning
 /** Resets the internal iteration state to the beginning of the current pool
  * without re-initializing or changing the pool itself. This allows iterating
  * through the same pool multiple times without the overhead of
  * re-initialization.
  *
  * This method should only be called if a pool has been initialized (via
  * init_random_pool() or init_representative_pool()). If no pool has been
  * initialized, this method should throw an exception.
  *
  * After calling reset_pool(), the next call to get_current_scenario() will
  * return the first scenario in the pool, and next_scenario() will iterate
  * from the beginning.
  *
  * @throws std::runtime_error If no pool has been initialized */

 virtual void reset_pool( void ) = 0;

/*--------------------------------------------------------------------------*/
 /// check if a pool has been initialized
 /** Returns true if either init_random_pool() or init_representative_pool()
  * has been called and the pool is ready for iteration. Returns false if no
  * pool has been initialized yet.
  *
  * This method is useful to check whether reset_pool() can be called instead
  * of initializing a new pool, which can be more efficient when the same pool
  * needs to be iterated multiple times.
  *
  * @return true if a pool is initialized and ready for iteration, false
  *         otherwise */

 [[nodiscard]] virtual bool is_pool_initialized( void ) const = 0;

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

 // type of the factory of ScenarioGenerator
 using ScenarioGeneratorFactory = boost::function< ScenarioGenerator *() >;

 // type of the map between strings and the factory of ScenarioGenerator
 using ScenarioGeneratorFactoryMap = std::map< std::string ,
                                               ScenarioGeneratorFactory >;

/*--------------------------------------------------------------------------*/
/*------------------------- PROTECTED METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Protected methods for handling static fields
 *
 * These methods allow derived classes to partake into static initialization
 * procedures performed once and for all at the start of the program. These
 * are typically related to factories.
 * @{ */

 /// method encapsulating the Solver factory
 /** This method returns the Solver factory, which is a static object. The
  * rationale for using a method is that this is the "Construct On First Use
  * Idiom" that solves the "static initialization order problem". */

 static ScenarioGeneratorFactoryMap & f_factory( void ) {
  static ScenarioGeneratorFactoryMap s_factory;
  return( s_factory );
  }

/*--------------------------------------------------------------------------*/
 /// empty placeholder for class-specific static initialization
 /** The method static_initialization() is an empty placeholder that is made
  * available to derived classes that need to perform some class-specific
  * static initialization besides these of any :ScenarioGenerator class,
  * i.e., the management of the factory. This method is invoked by the
  * SMSpp_insert_in_factory_cpp_* macros [see SMSTypedefs.h] during the
  * standard initialization procedures. If a derived class needs to perform
  * any static initialization, it just has to do this into its version of
  * this method; if not, it just has nothing to do, as the (empty) method of
  * the base class will be called.
  *
  * This mechanism has a potential drawback in that a redefined
  * static_initialization() may be called multiple times. Assume that a
  * derived class X redefines the method to perform something, and that a
  * further class Y is derived from X that has to do nothing, and that
  * therefore will not define Y::static_initialization(): them, within the
  * SMSpp_insert_in_factory_cpp_* of Y, X::static_initialization() will be
  * called again. This is not likely to happen to :ScenarioGenerator, and it
  * is easily managed in case. */

 static void static_initialization( void ) {}

/** @} ---------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
 // Definition of ScenarioGenerator::private_name() (pure virtual)

 [[nodiscard]] virtual const std::string & private_name( void ) const = 0;

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };   // end( class ScenarioGenerator )

/*--------------------------------------------------------------------------*/
/*----------------- CLASS MultiStageScenarioGenerator ----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// abstract base class for a handler of multistage scenarios
/** MultiStageScenarioGenerator derives from ScenarioGenerator and extends its
 * abstract interface from the single-stage one, where the stochastic is a
 * single random variable X_0 (although in fact this can be the composition of
 * different random variables, but all of them taking values at the same
 * instant in time, to the general framework of multistage problem where a
 * number T > 0 of time instant (stages) are defined and a scenario ( x_0 ,
 * ... , x_{T-1} ) is a realization of a finite horizon stochastic process
 * ( X_0 , ... , X_{T-1} ), where for each 0 <= t <= T - 1, the random
 * variable X_t has values in some Euclidean space R^(d_t). The base class
 * makes no assumption on the "internal structure" of x_t, which is just a
 * vector of the appropriate length d_t, whose internal structure is
 * supposed to be fixed and known to the user.
 *
 * In the following, at a stage 0 < t < T, we denote by H_t = ( X_0 = x0 ,
 * ... , X_{t-1} = x_{t-1} ) the history of the process up to time t
 * (H_0 being obviously empty). In general, in order to compute a scenario
 * ( x_0 , ... , x_{T-1} ) from the input stochastic process ( X_0 , ... ,
 * X_{T-1} ), the class will need to be able to draw for every t < T from the
 * random variable X_t | H_t. Yet, the class makes no specific assumption on
 * how the variable X_t depends on H_t. For instance, the case where each
 * stage is a random variable X_t independent of the past history H_t (i.e.,
 * X_t | H_t = X_t) is one case, where a sensible implementation is just as
 * a vector of ScenarioGenerator *, one for each of the stages. Other
 * possibilities exist, e.g., the Markovian setting where each stage depends
 * only from the results of the previous stage (that is, in fact, H_t =
 * ( X_{t-1} = x_{t-1} )) up to the general one where the random variables
 * depend on the entire history. The interfaces strive to avoid
 * distinguishing those cases by having each stage x_t generated / read-only
 * when its history H_t is clearly specified, irrespectively on what the
 * dependency is.
 *
 * Basically, the class provides means for moving along (back and forth) the
 * time dimension of the variable, since once the stage t is known in the
 * multistage setting then a scenario x_t for that stage "looks a lot like"
 * in the single-stage case, which is why a large part of the interface can
 * be inherited from the base ScenarioGenerator.
 *
 * See the comments on the base ScenarioGenerator class for more details. */

class MultiStageScenarioGenerator : public ScenarioGenerator
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
 *  @{ */

/// type for the index characterizing a specific size
/** Integer type that can be used to identify one of the possible stages t \in
 * 0, ..., T - 1, i.e., it has to be at least as large to hold T - 1. */

 using StageIndex = unsigned short;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// constexpr for "infinitely many stages"
 static constexpr StageIndex INFStage = Inf< StageIndex >();

/** @} ---------------------------------------------------------------------*/
/*------------ CONSTRUCTING AND DESTRUCTING ScenarioGenerator --------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing ScenarioGenerator
 *  @{ */

 /// constructor: initialize the few data structures of the base class

 MultiStageScenarioGenerator( void ) : ScenarioGenerator() {}

/*--------------------------------------------------------------------------*/
 /// destructor

 virtual ~MultiStageScenarioGenerator() = default;

/** @} ---------------------------------------------------------------------*/
/*- METHODS FOR READING THE STATIC DATA OF THE MultiStageScenarioGenerator -*/
/*--------------------------------------------------------------------------*/
/** @name Reading the static data of the MultiStageScenarioGenerator
 *
 * Note that in MultiStageScenarioGenerator get_scenario_size() does not
 * always report the same number, as what it is supposed to report is the
 * size of the scenario (the length of the vector) for the *current* time
 * instant.
 *  @{ */

 /// returns the number of stages

 [[nodiscard]] virtual StageIndex get_stage_number( void ) = 0;

/*--------------------------------------------------------------------------*/
 /// returns true if the stages are mutually independent
 /** A MultiStageScenarioGenerator is *stage-independent* when, for every t,
  * the random variable X_t does not depend on the history H_t. In this case
  * the scenario tree degenerates into a Cartesian product of single-stage
  * pools, one per stage, and the size of each per-stage pool can be set
  * independently. Consumers that exploit this structure (e.g.,
  * SDDPBlock::prepare_multi_stage_generator_pool() walking each stage on
  * its own) require this property; this method lets them check it
  * generically, without resorting to dynamic_cast.
  *
  * The base class returns false; only :MultiStageScenarioGenerator that
  * actually implement stage-independent scenarios should override and
  * return true. */

 [[nodiscard]] virtual bool is_stage_independent( void ) const {
  return( false );
  }

/** @} ---------------------------------------------------------------------*/
/*-------------------- METHODS FOR READING THE SCENARIOS -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Managing the time component while reading the scenarios
 *
 * MultiStageScenarioGenerator inherits the concepts of "pool" and "current
 * scenario" from ScenarioGenerator and adds that of "position in the
 * process". Because in general X_t depends on the history H_t, a position
 * is not just a stage index: it is the whole path ( x_0 , ... , x_{t-1} )
 * followed to reach it. A position is materialized by a View [see below],
 * which pins one history H and *is* a single-stage ScenarioGenerator over
 * the realizations that can be drawn next. Hence, all the "current
 * scenario" methods inherited from ScenarioGenerator (get_scenario_size(),
 * get_current_scenario(), next_scenario(), reset_pool(), init_*_pool(),
 * ...) are meant to be called on a View, and refer to the position that
 * View pins. The MultiStageScenarioGenerator is a ScenarioGenerator itself,
 * and as such it behaves exactly like the View pinned at the root, i.e., it
 * reads the realizations of the first random variable; this is a
 * convenience for the frequent case where only those are needed.
 *
 * The View is the *only* mechanism to move along the time dimension:
 * View::descend() fixes the realization currently selected and moves to the
 * next stage, View::climb() moves back to the previous one. It is important
 * to stress that drawing another realization of X_t with the same history
 * H_t is *not* descend() but the inherited next_scenario(): descend() moves
 * to a new *partial* scenario ( H_t , x_t ), i.e., it extends the history.
 * One only gets a "full" scenario ( x_0 , ... , x_{T-1} ) after having
 * called get_current_scenario() at each position of a root-to-leaf path.
 * Hence, when x_t is drawn from X_t the history H_t, if necessary, is
 * known.
 *
 * A consequence of this, however, is that the parameter size in init_*_pool()
 * is only a "rough" description of the scenario set. The bound is on the
 * total number of "full" scenarios, and it does not say how many of those
 * share a common history and how. For instance, if the random variables X_t
 * are independent, one can easily obtain a set of scenarios by deciding the
 * number s_t of different realizations for each stage t, and then size
 * = s_0 * ... * s_{T-1}; however, clearly there can be many different ways
 * to choose the s_t that give the same size. Thus, in this case one could
 * have get_support_size() to return s_t for the stage of the View,
 * similarly to get_scenario_size(), and even specify the size of each
 * sub-pool. However, there does not seem to be a general way of doing this
 * for all reasonable distributions that actually depend on H_t, in that the
 * right number of branching in the scenario tree at a node may depend on the
 * whole path to the root. The design choice is therefore that these more
 * detailed information about how the scenario tree is organized depends on
 * the specific derived class of MultiStageScenarioGenerator and read together
 * with its description in deserialize(); hence, they cannot be changed on the
 * fly. If this turns out to be a problem, then some Configuration will have
 * to be used, but we'll cross that bridge when (if) we come to it.
 *
 * Anyway, this fully justifies the return value of next_scenario(), since in
 * MultiStageScenarioGenerator, as opposed to its base class, it is not known
 * beforehand how many sub-scenarios sharing the same history H_t will be
 * available. */

 // Semantics of the inherited pool methods on
 // :MultiStageScenarioGenerator. Because the base ScenarioGenerator has
 // no concept of stage, the refinement of its API can only be documented
 // here.
 //
 // The inherited scalar init_*_pool( K ) overloads from ScenarioGenerator
 // act on the pool of the position they are called at, mirroring the
 // single-stage semantics of get_scenario_size() /
 // get_current_scenario() / next_scenario(). To (re-)init the pool of
 // more than one stage, the caller walks the tree with View::descend()
 // and issues a fresh init_*_pool( K_t ) at each. No per-stage or vector
 // overloads are provided: the explicit walk is the canonical pattern and
 // keeps the API of ScenarioGenerator and MultiStageScenarioGenerator
 // identical.
 //
 // Symmetrically, the inherited reset_pool() rewinds
 // get_current_scenario() to the first realization available at the
 // position it is called at, and does not move that position: rewinding
 // the position itself is done with View::climb(), or by taking a fresh
 // root_view().

/** @} ---------------------------------------------------------------------*/
/*------------------- VIEW-BASED (HISTORY-PINNED) ACCESS -------------------*/
/*--------------------------------------------------------------------------*/
/** @name View-based, history-pinned access
 *  @{ */

 /// a lightweight, history-pinned position in the process, itself a SG
 /** A View pins a history H (a sub-path from the root down to a node) and
  * *is* a single-stage ScenarioGenerator over the realizations that can be
  * drawn next, i.e. over X_{t+1} | H, with their conditional probabilities
  * as the pool weights. A consumer that only needs "the scenarios at this
  * point" (e.g. a TwoStageStochasticBlock) thus uses a View as an ordinary
  * ScenarioGenerator, unaware of sitting inside a larger tree.
  *
  * A View is the only mechanism to move along the time dimension: descend()
  * fixes the realization currently selected and moves one stage forward,
  * climb() moves back to the previous one, and each undoes the other, so
  * that a single View freely roams the whole tree. Both move *this* View: to
  * retain a position while moving on, take a clone() of it first.
  *
  * A View is meant to be a cheap handle (conceptually a pointer to the
  * generator plus the pinned position), so that arbitrarily many can be held
  * at once, e.g. one per node to build the corresponding Block tree.
  * Whether distinct Views can also be *used* concurrently depends on the
  * specific derived class: it is the case whenever the underlying data are
  * immutable and each View carries its own cursor, but a generator whose
  * positions share a mutable cursor cannot promise it.
  *
  * Every method of a View but the init_*_pool() ones only *reads* the
  * process: get_support_size(), get_current_scenario(), next_scenario(),
  * reset_pool(), descend(), climb() and clone() leave the tree exactly as
  * they found it, and can therefore be used by arbitrarily many Views at
  * once. The init_*_pool() ones are the only ones that *write*: called on a
  * View they (re-)define the pool of the node it is pinned at, which in a
  * generator that constructs the tree on the fly is where the children of
  * that node come into existence. The pool of a View must therefore have
  * been initialised [see is_pool_initialized()] before it can be read or
  * descend()-ed into, while a generator whose tree is entirely materialised
  * upfront simply has all its Views born initialised.
  *
  * Since init_*_pool() re-defines a pool, it invalidates the Views that
  * were relying on it: those pinned at the very same node, whose pool has
  * just changed under their feet, and, unless is_stage_independent() is
  * true, those pinned anywhere below it, whose sub-tree has just been
  * regenerated. Views pinned elsewhere are unaffected: those on a different
  * sub-tree share no node with it, and those above it only see their own
  * children, which have not changed. Note that this depends only on where
  * the Views are pinned now, not on the order in which they were created,
  * since climb() allows a View to be anywhere at any time. A View knows
  * whether it has been invalidated [see is_valid()] and refuses to be used
  * if it has, rather than silently reading stale data; the View that
  * init_*_pool() was called on obviously remains valid. */

 class View : public ScenarioGenerator
 {
  public:

   /// move into the realization currently selected, one stage forward
   /** Moves this View one stage forward, extending the pinned history H with
    * the realization currently selected by the inherited
    * get_current_scenario(); its pool becomes that of X_{t+1} | ( H , x_t ),
    * with the first realization selected. Returns false, leaving the View
    * unchanged, if there is no further stage, i.e., the realization
    * currently selected is a leaf of the scenario tree.
    *
    * Returning false is reserved for that legitimate answer: descending
    * from a View whose pool has not been initialised yet [see
    * is_pool_initialized()] is instead a usage error, since which the
    * children are is precisely what a pool defines, and throws exception.
    * The arrival View is initialised, or not, exactly as any other View
    * pinned at that node would be. */

   virtual bool descend( void ) = 0;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
   /// move back to the position this View has descended from
   /** Moves this View one stage backwards, dropping the last realization of
    * the pinned history H; its pool becomes that of the previous stage, with
    * the realization that had been descended into selected, so that climb()
    * undoes descend() and vice-versa. Returns false, leaving the View
    * unchanged, if it is pinned at the root, where H is empty. */

   virtual bool climb( void ) = 0;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
   /// returns an independent copy of this View, pinned at the same position
   /** Returns a View that pins the same history H as this one and selects the
    * same realization in the pool, but has its own cursor: the two then move
    * independently. This is how one retains a position while moving on, and
    * how a Views is handed out to each of several concurrent consumers. */

   [[nodiscard]] virtual std::unique_ptr< View > clone( void ) const = 0;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
   /// the stage (depth of the pinned history H) of this View

   [[nodiscard]] virtual StageIndex stage( void ) const = 0;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
   /// true if this View still refers to the process it was created on
   /** Returns true if this View can still be used, and false if some
    * init_*_pool() has meanwhile re-defined a pool it was relying on [see
    * the general notes of this class for which Views a pool re-definition
    * invalidates]. An invalidated View is inert: every other method of it
    * throws exception rather than reading data that is no longer there, so
    * that an owner that keeps Views around across pool re-definitions can
    * ask, and rebuild the ones it has lost, instead of finding out much
    * later. A View of a generator that never re-defines a pool is of course
    * valid for its whole life. */

   [[nodiscard]] virtual bool is_valid( void ) const = 0;

 };   // end( class View )

/*--------------------------------------------------------------------------*/
 /// create a View pinned at the root (empty history, stage 0)
 /** Returns a View pinned at the root of the scenario tree: its history H is
  * empty and its pool are the realizations of the first random variable.
  * From it the whole tree is reached via View::descend() and View::climb().
  * Since the View is the only mechanism to read a scenario tree, this is the
  * entry point every :MultiStageScenarioGenerator has to provide. */

 [[nodiscard]] virtual std::unique_ptr< View > root_view( void ) const = 0;

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED FIELDS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

 // probably not needed SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };   // end( class MultiStageScenarioGenerator )

/** @} end( group( Classes in ScenarioGenerator ) ) */
/*--------------------------------------------------------------------------*/
/*------------------- inline methods implementation ------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* ScenarioGenerator.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File ScenarioGenerator.h ------------------------*/
/*--------------------------------------------------------------------------*/

