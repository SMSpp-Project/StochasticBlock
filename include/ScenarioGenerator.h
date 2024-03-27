/*--------------------------------------------------------------------------*/
/*--------------------- File ScenarioGenerator.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the definition of the *abstract* classes ScenarioGenerator
 * and MultiStageScenarioGenerator, that represent the interface between a user
 * that needs to generate instantiations of random processes (scenarios),
 * typically for the purpose of forming / solving stochastic optimization
 * problems, and the actual data structures holding the data that contains /
 * generates these instantiations.
 *
 * There can be several different ways in which such data is organized, starting
 * from the fact that there is only one "stage" (all the random variables are
 * revealed at the same time) or multiple ones (there is a different random
 * variable for each of a discrete set of points in time, called stages, and
 * random variables are progressively revealed as time passes through them).
 * Once the point in time is fixed, though, the data looks the same in both
 * cases; this is why the base class ScenarioGenerator defines the interface for
 * retrieving the data concerning a single stage, and
 * MultiStageScenarioGenerator derives from it by adding the concepts related to
 * handling the different points in time (stages).
 *
 * Apart from that, it is expected that different concrete derived classes will
 * be implemented to consider relevant special cases, such as:
 *
 * - DiscreteScenarioSet deriving from ScenarioGenerator for the case where the
 *   (single-stage) scenarios are just the result of a simple unique random
 *   variable with a discrete distribution (one element for each scenario, just
 *   provided in input in their final form).
 *
 * - MultiStageDiscreteScenarioSet deriving from MultiStageScenarioGenerator for
 *   the case where the (multistage) scenarios are just the result of a simple
 *   unique random variable with a discrete distribution (one element for each
 *   scenario, just provided in input in their final form).
 *
 * - IndependentMultiStageScenarioGenerator where each stage is a random
 *   variable independent from the past history; a possible implementation is
 *   just as a vector of ScenarioGenerator *, one for each of the stages.
 *
 * Many other possibilities exist, e.g., single-stage scenarios constructed from
 * the combination of an arbitrary number of distributions (discrete, continuous
 * of various types, ...), multi-stage scenarios in the Markovian setting where
 * each stage is (an arbitrary combination of random variables) depending only
 * from the results of the previous stage, general multi-stage scenarios where
 * the random variables depend on the entire previous history, and so on. In
 * many cases the user is not concerned with these details, and therefore they
 * can just program agaist the [MultiStage]ScenarioGenerator interface.
 *
 * \author Antonio Frangioni \n Dipartimento di Informatica \n Universita' di
 *         Pisa \n
 *
 * \author Benoit Tran \n Dipartimento di Informatica \n Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
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
#include <span> // To be added in SMSTypedefs.h includes of the std library?


/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it {

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
 * address in part both. This is because both share the interface that allows to
 * access the data (realization) of one specific point in time (stage), of which
 * there is only one in the case of ScenarioGenerator.
 *
 * In the general framework of a multistage problem, a number T > 0 of time
 * instant (stages) are defined. Then, a scenario ( x_0 , ... , x_{T-1} ) is a
 * realization of a finite horizon stochastic process ( X_0 , ... , X_{T-1} ),
 * where for each 0 <= t <= T - 1, the random variable X_t has values in some
 * euclidean space R^(d_t). The base class makes no assumption on the "internal
 * structure" of x_t, which can in general be the composition of different
 * random variables (say, a pair ( c , a ) where c is the vector c[ j ] denoting
 * the cost of a given commodity at a predetermined set of locations j \in J,
 * and a is the vector a[ j ] denoting the availability of the same commodity at
 * the same places; or the matrix of pairs where each row correspond to a
 * different commodity, ...). As far as the base classes are concerned, a
 * scenario x_t for a given stage t is just a vector of appropriate length d_t,
 * whose internal structure is supposed to be fixed and known to the user. Also,
 * the base classes make no assumption on how the variable X_t depends on ( x_0
 * , ... , x_{t-1} ) (the history at stage t).
 *
 * When T = 1, the time dimension is implicit, i.e., the stochastic process is
 * its only element, the random variable X_0, and a scenario is its unique value
 * x_0. This is why there are two separate classes: ScenarioGenerator and
 * MultiStageScenarioGenerator. In the former the time dimension is fixed, and
 * therefore there is no need of methods to manage it, differently from the
 * latter. However, once the stage t is known in the multistage setting then a
 * scenario x_t for that stage "looks a lot like" in the single-stage case,
 * which is why a large part of the interface is shared.
 *
 * Other notable elements of ScenarioGenerator:
 *
 * - it has a factory (like, Block, Solver, Configuration, ...) managed by means
 *   of the standard SMSpp_insert_in_factory_h and
 *   SMSpp_insert_in_factory_cpp_[0_t]() macros;
 *
 * - it therefore has a static method new_ScenarioGenerator( std::string ) that
 *   constructs a :ScenarioGenerator of the specific type dictated by the string
 *   argument (using the factory)
 *
 * - it has an in-built hierarchy of deserialization methods, starting with
 *   deserialize( std::string ) that deserializes a :ScenarioGenerator from the
 *   netCDF file with the given filename, passing through new_ScenarioGenerator(
 *   netCDF::NcGroup ) that reads some top-level information from the group
 *   (either the class type or a filename), finally down to the pure virtual
 *   deserialize( netCDF::NcGroup ) that actually implements the deserialization
 *   for the specific :ScenarioGenerator
 *
 * - it allows to set the seed of the pseudo-random number generator(s) possibly
 *   involved in the (sampling) operations
 *
 * - distinguishes different types of "scenario pools" supporting different use
 *   cases
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
/** Type for reading the data characterizing a scenario. The base class does not
 * make any assumption on any "structure" that a scenario may have; a scenario
 * is just a vector of numbers (double). It is assumed that the representation
 * of the "current scenario" that can be read as a contiguous array of double
 * will always be available inside any implementation of the abstract base
 * class, hence a Scenario is just a std::span<> provided read-only access to
 * that data structure. */
 
 using Scenario = std::span< const double >;

/*--------------------------------------------------------------------------*/
/// type for the index characterizing a specific scenario
/** Integer type that can be used to identify one of the possible finite number
 * of realizations that can be produced. Basically implies the maximum number of
 * them. */

 using ScenarioIndex = unsigned int;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// constexpr for "infinitely many scenarios"
 static constexpr ScenarioIndex INFScenario = Inf< ScenarioIndex >(); 

/*--------------------------------------------------------------------------*/
/// type for the scenario size (the length of the double vector)
/** Integer type that can be used to specify the length of the double vector
 * containing the data of the scenario. */

 using ScenarioSize = size_t;

/** @} ---------------------------------------------------------------------*/
/*------------ CONSTRUCTING AND DESTRUCTING ScenarioGenerator --------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing ScenarioGenerator
 *  @{ */

 /// constructor: initialise the few data structures of the base class
 ScenarioGenerator( void ) { }

/*--------------------------------------------------------------------------*/
 /// construct a :ScenarioGenerator of specific type using the factory
 /** Use the ScenarioGenerator factory to construct a :ScenarioGenerator object
  * of type specified by classname (a std::string with the name of the class
  * inside). Note that the method is static because the factory is static, hence
  * it is to be called as
  *
  *   ScenarioGenerator *mySG = ScenarioGenerator::new_ScenarioGenerator(
  *                   some_class );
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
  * @param classname The name of the :Solver class that must be constructed. */

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
  * second-level method to allow indirect file access, see the comment there for
  * the format of the file.
  *
  * If anything goes wrong with deserialization, nullptr is returned.
  *
  * Note that the method is static, hence it is to be called as
  *
  *      auto mySG = ScenarioGenerator::deserialize( some_file );
  *
  * i.e., without any reference to any specific ScenarioGenerator. */

[[nodiscard]] static ScenarioGenerator* deserialize(const std::string& filename) {
    try {
      // Relies on netCDF API to handle is_open() check
      netCDF::NcFile dataFile(filename, netCDF::NcFile::read);

      return( ScenarioGenerator::new_ScenarioGenerator( dataFile ) );
    }
    catch (netCDF::exceptions::NcException& e) {
        std::cerr << "Error opening or processing netCDF file: " << e.what() << std::endl;
    }
    return nullptr;
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
  * - A "direct" group that contains at least the string attribute "type"; this
  *   is used it in the factory to construct an "empty" :ScenarioGenerator of
  *   that type [see new_ScenarioGenerator( std::string & )], and then the
  *   method deserialize( netCDF::NcGroup ) of the newly minted
  *   :ScenarioGenerator is invoked (with argument \p group) to finish the work.
  *
  * - An "indirect" group that just need to contain the single string attribute
  *   "filename"; in this case, the attribute is used as argument for a call to
  *   deserialize( const std::string & ) that will extract the
  *   :ScenarioGenerator by the corresponding netCDF file.
  *
  * In case \p group contains both "type" and "filename", the first takes the
  * precedence (direct groups have precedence over indirect ones).
  *
  * Note that this method is static (see the previous versions for comments
  * about it) and returns a pointer to ScenarioGenerator, hence it has to have a
  * different name from deserialize( netCDF::NcGroup ) (since the signature is
  * the same but for the return type).
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
  * A group containing a :Block must have a mandatory string attribute "type"
  * that contains the classname() of the :ScenarioGenerator, which is actually
  * useful at the higher level of the deserialize() hierarchy where the
  * :ScenarioGenerator has already to be constructed, rather than at this point
  * where it clearly already has, plus whatever other information is required by
  * the specific :ScenarioGenerator.
  *
  *      THIS IS THE METHOD TO BE IMPLEMENTED BY DERIVED CLASSES
  *
  * and in fact it is pure virtual. */
 virtual void deserialize( const netCDF::NcGroup & group ) = 0;

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
  * can ever output. This is done so that the user avoids to call the
  * init_*_pool() methods [see] with requiring too high a number of scenarios.
  * In the case of the random variable being continuous the support size is
  * (theoretically) infinite, so INFScenario is reported. This is done by the
  * base class implementation, so that the method only needs to be redefined by
  * derived classes implementing variables with finite support.
  */

 [[nodiscard]] virtual ScenarioIndex get_support_size( void ) {
  return( INFScenario );
  }

/*--------------------------------------------------------------------------*/
 /// returns the size of any scenario (length of the double vector)
 /** This method returns size of the double vector containing the data that
  * specify the instantiation x_0 of the random variable X_0. This will be the
  * size of the std::span< const double > (Scenario) returned by
  * get_current_scenario(), */
  
 [[nodiscard]] virtual ScenarioSize get_scenario_size( void ) = 0;

/*--------------------------------------------------------------------------*/
 /// getting the classname of this ScenarioGenerator
 /** Given a :ScenarioGenerator (i.e., any class derived by it), this method
  * returns a string with its class name; unlike std::type_info.name(), there
  * *are* guarantees, i.e., the name will always be the same.
  *
  * The method works by dispatching the private virtual method private_name().
  * The latter is automatically implemented by the SMSpp_insert_in_factory_cpp_*
  * macros [see SMSTypedefs.h], hence this comes at no cost since these have to
  * be called somewhere to ensure that any :ScenarioGenerator will be added to
  * the factory. Actually, since ScenarioGenerator::private_name() is pure
  * virtual, this ensures that it is not possible to forget to call
  * SMSpp_insert_in_factory_cpp_* for any :ScenarioGenerator because otherwise
  * it is a pure virtual class (unless the programmer purposely defines
  * private_name() without calling the macro, which seems rather pointless). */

 [[nodiscard]] const std::string & classname( void ) const {
  return( private_name() );
  }

/** @} ---------------------------------------------------------------------*/
/*-------------------- METHODS FOR READING THE SCENARIOS -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Generating the scenario pool and reading the scenarios
 *
 * ScenarioGenerator works based on the concept of "scenario pool". The user has
 * to request the (logical) creation of a finite set of scenarios that they can
 * then read. Scenario pools can be constructed differently with different aims
 * in mind (and a different cost). Once a scenario pool is initialised it can be
 * traversed by going through the scenarios in the pool one by one; at each
 * stage a current scenario is defined that can be read. Once the current
 * scenario is moved forward the previous current scenario is in principle no
 * longer available.
 *
 * The user is supposedly incapable of differentiating the scenarios of the pool
 * from one another, which means that there is no assumption on the order in
 * which the scenarios of the pool is visited. Also, there should be no need to
 * re-visit a scenario once the user has been visited it and moved on, so "the
 * clock of the pool only moves forward" (if the user needs to look at more than
 * one scenario at the time they will have to save them in their own data
 * structure.
 *  @{ */

 /// set the seed of the internal random generator
 /** Creating pools can well involve random sampling operations, and therefore
  * the use of pseudo-random numbers. To guarantee reproducibility, the seed for
  * the pseudo-random generator can be set with this method. It is conceivable
  * that some "complicated" ScenarioGenerator may require more than one
  * generator, if this is the case the method may have to be changed (but it is
  * not impossible to use the same seed for all of them, or to use the seed
  * within a random number generator that is used to generate the seeds, so just
  * one seed may well be enough). */

 virtual void set_seed( unsigned long seed = 0 ) = 0;

/*--------------------------------------------------------------------------*/
 /// generate a random pool of given size
 /** This method (logically) constructs the scenario pool from realizations of
  * the random variable, typically by sampling. One would expect the scenarios
  * to be uniformly distributed across the (possibly, infinite) set of possible
  * realizations, so that the corresponding scenarios can be used, e.g., during
  * a simulation phase to evaluate the quality of the decisions taken by some
  * optimization model that has only "seen" a (much) smaller number of scenarios
  * [see init_representative_pool()].
  *
  * \p size is the number of scenarios in the produced pool.
  *
  * Once the method returns, the first scenario in the pool is available as the
  * current scenario and can be read [see get_current_scenario()] and
  * get_current_scenario_probability()] right away. */

 [[nodiscard]] virtual void init_random_pool( ScenarioIndex size ) = 0;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate a "representative" pool of given size
 /** This method (logically) constructs the scenario pool from realizations of
  * the random variable, working hard to ensure that the set is "as much
  * representative as possible" of the whole (possibly, infinite) set of
  * possible realizations, so that the corresponding scenarios can be used,
  * e.g., to build an optimization model of "reasonable size" (so that it can be
  * solved) that still "sees enough of the random variable to take good
  * decisions". This may be nontrivial, in particular if the underlying variable
  * is discrete in nature and therefore something like an Optimal Transport
  * Problem need be solved. Thus, this method may have a nontrivial
  * computational cost. Then, the decisions taken by such model can be evaluated
  * with a simulation on a (much) larger number of "uniformly distributed"
  * scenarios [see init_random_pool()].
  *
  * \p size is the number of scenarios in the produced pool.
  *
  * Once the method returns, the first scenario in the pool is available as the
  * current scenario and can be read [see get_current_scenario()] and
  * get_current_scenario_probability()] right away. */

 [[nodiscard]] virtual void init_representative_pool(
						    ScenarioIndex size ) = 0;

/*--------------------------------------------------------------------------*/
 /// read the data of the current scenario
 /** Gets a std::span< const double > (Scenario) of size get_scenario_size()
  * [see] that contains the data that specify the instantiation x_0 of the
  * random variable X_0. This points to an internal data structure of the
  * :ScenarioGenerator, which is guaranteed to remain available up until the
  * next call of either next_scenario() or init_*_pool(), after which the data
  * structure may no longer be available and the span should no longer be used.
  * */

 [[nodiscard]] virtual Scenario get_current_scenario( void ) = 0;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// return the probability associated to the current scenario
 /** Does it always make sense??
  */
 
 [[nodiscard]] virtual double get_current_scenario_probability( void )  = 0;

/*--------------------------------------------------------------------------*/
 /// move the current scenario to the next scenario in the pool
 /** Tries to move forward the current scenario to the next scenario in the
  * pool, returns true if this succeeds, which should always happen if the
  * number of scenarios in the pool seen so far is less than the size specified
  * in init_*_pool(). Upon receiving a true response the user can call again
  * get_current_scenario[_probability]() to retrieve the data of the ew
  * scenario, knowing that any Scenario corresponding to the previous scenarios
  * seen has potentially been invalidated and should no longer be used. */

 [[nodiscard]] virtual bool next_scenario( void ) = 0;

/** @} ---------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

 using ScenarioGeneratorFactory = boost::function< ScenarioGenerator *() >;
 // type of the factory of ScenarioGenerator

 using ScenarioGeneratorFactoryMap = std::map< std::string ,
                                               ScenarioGeneratorFactory >;
 // type of the map between strings and the factory of ScenarioGenerator

/*--------------------------------------------------------------------------*/
/*------------------------- PROTECTED METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Protected methods for handling static fields
 *
 * These methods allow derived classes to partake into static initialization
 * procedures performed once and for all at the start of the program. These are
 * typically related with factories.
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
 /** The method static_initialization() is an empty placeholder which is made
  * available to derived classes that need to perform some class-specific static
  * initialization besides these of any :ScenarioGenerator class, i.e., the
  * management of the factory. This method is invoked by the
  * SMSpp_insert_in_factory_cpp_* macros [see SMSTypedefs.h] during the standard
  * initialization procedures. If a derived class needs to perform any static
  * initialization it just have to do this into its version of this method; if
  * not it just has nothing to do, as the (empty) method of the base class will
  * be called.
  *
  * This mechanism has a potential drawback in that a redefined
  * static_initialization() may be called multiple times. Assume that a derived
  * class X redefines the method to perform something, and that a further class
  * Y is derived from X that has to do nothing, and that therefore will not
  * define Y::static_initialization(): them, within the
  * SMSpp_insert_in_factory_cpp_* of Y, X::static_initialization() will be
  * called again. This is not likely to happen to :ScenarioGenerator and it is
  * easily managed in case. */

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
 * different random variables, but all of them taking values at the same instant
 * in time, to the the general framework of multistage problem where a number T
 * > 0 of time instant (stages) are defined and a scenario ( x_0 , ... , x_{T-1}
 * ) is a realization of a finite horizon stochastic process ( X_0 , ... ,
 * X_{T-1} ), where for each 0 <= t <= T - 1, the random variable X_t has values
 * in some euclidean space R^(d_t). The base class makes no assumption on the
 * "internal structure" of x_t, which is just a vector of appropriate length
 * d_t, whose internal structure is supposed to be fixed and known to the user.
 *
 * In the following, at a stage 0 < t < T, we denote by H_t = ( X_0 = x0 , ... ,
 * X_{t-1} = x_{t-1} ) the history of the process up to time t (H_0 being
 * obviously empty). In general, in order to compute a scenario ( x_0 , ... ,
 * x_{T-1} ) from the input stochastic process ( X_0 , ... , X_{T-1} ), the
 * class will need to be able to draw for every t < T from the random variable
 * X_t | H_t. Yet, the class makes no specific assumption on how the variable
 * X_t depends on H_t. For instance, the case where each stage is a random
 * variable X_t independent from the past history H_t (i.e., X_t | H_t = X_t) is
 * one case, where a sensible implementation is just as a vector of
 * ScenarioGenerator *, one for each of the stages. Other possibilities exist,
 * e.g., the Markovian setting where each stage depends only from the results of
 * the previous stage (that is, in fact H_t = ( X_{t-1} = x_{t-1} )) up to the
 * general one where the random variables depend on the entire history. The
 * interfaces strives to avoid to distinguish those cases by having each stage
 * x_t generated / read only when its history H_t is clearly specified,
 * irrespectively on what the dependency is.
 *
 * Basically, the class provide means for moving along (back and forth) the time
 * dimension of the variable, since once the stage t is known in the multistage
 * setting then a scenario x_t for that stage "looks a lot like" in the
 * single-stage case, which is why a large part of the interface can be
 * inherited from the base ScenarioGenerator.
 *
 * See the comments to the base ScenarioGenerator class for more details. */

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

 /// constructor: initialise the few data structures of the base class

 MultiStageScenarioGenerator( void ) : ScenarioGenerator() {}

/*--------------------------------------------------------------------------*/
 /// destructor

 virtual ~MultiStageScenarioGenerator() = default; 

/** @} ---------------------------------------------------------------------*/
/*- METHODS FOR READING THE STATIC DATA OF THE MultiStageScenarioGenerator -*/
/*--------------------------------------------------------------------------*/
/** @name Reading the static data of the MultiStageScenarioGenerator
 *
 * Note that in MultiStageScenarioGenerator get_scenario_size() does not always
 * report the same number, as what it is supposed to report is the size of the
 * scenario (the length of the vector) for the *current* time instant.
 *  @{ */

 /// returns the number of stages

 [[nodiscard]] virtual StageIndex get_stage_number( void ) = 0;

/** @} ---------------------------------------------------------------------*/
/*-------------------- METHODS FOR READING THE SCENARIOS -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Managing the time component while reading the scenarios
 *
 * MultiStageScenarioGenerator inherits the concepts of "pool" and "current
 * scenario" from ScenarioGenerator and adds that of "current stage". The
 * current stage is initialised to t = 0 (first stage, with empty H_t) when the
 * current scenario is (logically) created, i.e., right after the call to
 * init_*_pool() or that of next_scenario(), and either increased with
 * next_stage() or decreased by previous_stage(). Hence, when
 * get_current_scenario() is called, it refers to the scenario *for the current
 * stage*, and so does get_scenario_size() (although the same information can be
 * inferred from the returned Scenario.
 *
 * It is important to stress that, in MultiStageScenarioGenerator, the meaning
 * of next_stage() is to draw (if possible) another realization of X_t, with t =
 * get_current_stage(). Hence, this creates a new *partial* scenario ( H_t , x_t
 * ) with the same H_t as before the call and a different x_t. One only gets a
 * "full" scenario ( x_0 , ... , x_{T-1} ) after having called
 * get_current_scenario() for all t = 0, ..., T - 1. Hence, when x_t is drawn
 * from X_t the history H_t, if necessary, is known.
 *
 * A consequence of this, however, is that the parameter size in init_*_pool()
 * is only a "rough" description of the scenario set. The bound is on the total
 * number of "full" scenarios, and it does not say how many of those share a
 * common history and how. For instance, if the random variables X_t are
 * independent, one can easily obtain a set of scenarios by deciding the number
 * s_t of different realizations for each stage t, and then size = s_0 * ... *
 * s_{T-1}; however, clearly there can be many different ways to choose the s_t
 * that give the same size. Thus, in this case one could have get_support_size()
 * to return s_t for t = get_current_stage(), similarly to get_scenario_size(),
 * and even specify the size of each sub-pool. However, there does not seem to
 * ba a general way of doing this for all reasonable distributions that actually
 * depend on H_t, in that the right number of branching in the scenario tree at
 * a node may depend on the whole path to the root. The design choice is
 * therefore that these more detailed information about how the scenario tree is
 * organized are depending on the specific derived class of
 * MultiStageScenarioGenerator and read together with its description in
 * deserialize(); hence, they cannot be changed on the fly. If this turns out to
 * be a problem than some Configuration will have to be used, but we'll cross
 * that bridge when (if) we'll come to it.
 *
 * Anyway, this fully justifies the return value of next_scenario(), since in
 * MultiStageScenarioGenerator, as opposed to its base class, it is not known
 * beforehand how many sub-scenarios sharing the same history H_t will be
 * available. */

 /// get the current stage t, a number in 0, ..., get_stage_number() - 1

 [[nodiscard]] virtual StageIndex get_current_stage( void ) = 0;

/*--------------------------------------------------------------------------*/
 /// increase the current stage by one (if possible)
 /** Increase the value of the current stage t = get_current_stage() to t' =
  * min( t + 1 , T - 1 ); returns true if it succeeds, i.e., t < t'. If true is
  * returned, then this "resets the current scenario", i.e., the Scenario
  * (std::range< const double >) returned by the latest call to
  * get_current_scenario() is invalidated and the first scenario x_{t'} of the
  * next random variable X_{t'} is drawn, ideally subject the current history (
  * x_0 , ... , x_t = x_{t' - 1} ), and made available via the next call to
  * get_current_scenario(). */
  
 [[nodiscard]] virtual bool next_stage( void ) = 0;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// decrease the current stage by \p step (if possible)
 /** Decrease the value of the current stage t = get_current_stage() to t' =
  * max( t - step , 0 ); returns true if it succeeds, i.e., t' < t. If true is
  * returned, then this "resets the current scenario", i.e., the Scenario
  * (std::range< const double >) returned by the latest call to
  * get_current_scenario() is invalidated and the first scenario x_{t'} of the
  * next random variable X_{t'} is drawn, ideally subject the current history (
  * x_0 , ... , x_{t' - 1} ), if nonempty, and made available via the next call
  * to get_current_scenario(). */

 virtual void previous_stage( StageIndex step = 1 ) = 0;

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

