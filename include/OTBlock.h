/*--------------------------------------------------------------------------*/
/*-------------------------- File OTBlock.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class OTBlock, which implements the Block
 * concept [see Block.h] for the Optimal Transport problem.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 * 
 * \author Romain Pujol \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 * 
 * \author Benoît Tran \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Romain Pujol and Benoît Tran
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __OTBlock
 #define __OTBlock  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"

#include "LinearFunction.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "OneVarConstraint.h"

#include "Solution.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
 class OTBlock;     // forward declaration of OTBlock

 class OTSolution;  // forward declaration of OTSolution

/*--------------------------------------------------------------------------*/
/*----------------------- OTBlock-RELATED TYPES ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup OTBlock_TYPES OTBlock-related types
 *  @{ */

 using p_OTBlock = OTBlock *;  ///< a pointer to OTBlock

 using Vec_OTBlock = std::vector< p_OTBlock >;
 ///< a vector of pointers to OTBlock

 using Vec_OTBlock_it = Vec_OTBlock::iterator;
 ///< iterator for a Vec_OTBlock

 using c_Vec_OTBlock = const Vec_OTBlock;
 ///< a const vector of pointers to OTBlock

 using c_Vec_OTBlock_it = c_Vec_OTBlock::iterator;
 ///< iterator for a c_Vec_OTBlock

/** @}  end( group( OTBlock_TYPES ) ) */ 
/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup OTBlock_CLASSES Classes in OTBlock.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS OTBlock --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for the Optimal Transport problem
/** The OTBlock class implements the Block concept [see Block.h] for the
 * Optial Transport (OT) problem.
 * 
 * The data of the problem is two discrete probability distributions P and Q. 
 * Continuous distributions are not handled by OTBlock. The probability
 * distributions P = (p_i)_{i \in N} and Q = (q_i)_{i \in M} are assumed to 
 * respectively have n = |N| atoms (x_i)_{i \in N} and m = |M| atoms 
 * (y_i)_{i \in M}. The spaces in which the atoms of P and Q belong are not 
 * important but OT requires the knowledge of a cost function that associates a 
 * positive real c_{ij} to each atom i of P to each atom j of Q which stands
 * for the cost of moving the atom i to the atom j. The real c_{ij} is also 
 * the cost of moving atom j to atom i.
 * 
 * Informally, the objective of OT is to find the best transport plan 
 * \pi = (\pi_{ij})_{ij} such that the sum of the costs of moving the atoms of P
 * to the atoms of Q is minimized. Formally, the optimization problem we aim to 
 * solve is:
 * \f[
 *  \min_\pi \sum_{ (i, j) \in N x M } c[ i , j ] \pi[ i, j ]
 * \f]
 * \f[
 *  \sum_{ j \in M } \pi[ i , j ] = p_i \quad i \in N             (1)
 *  \sum_{ i \in N } \pi[ i, j ] = q_i \quad j \in M              (2)
 *  \pi[ i , j ] \geq 0.                                          (3)
 * \f]
 * 
 * The n+m equations in (1) and (2) are the mass conservation constraints. 
 * Along with the n*m equations in (3), the constraints of the OT problem 
 * encodes the fact that the transport plan \pi is made of 
 * weights \pi_{ij} that represents the proportion of the i-th atom of P sent to
 * the j-th atom of Q and vice versa. 
 * 
 * TODO: ADD A REMARK ON STATIC AND DYNAMIC VARIABLES WHEN IT IS CLEARER POST 
 * DISCUSSION WITH ROMAIN
 *
*/

class OTBlock : public Block
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * TODO:
 * 
 @{ */

/*--------------------------------------------------------------------------*/

 typedef double TNumber;                     ///< type of transport plan component
 typedef const TNumber c_TNumber;            ///< a read-only TNumber

 typedef std::vector< TNumber > Vec_TNumber; ///< a vector of TNumber
 typedef const Vec_TNumber c_Vec_TNumber;    ///< a const vector of TNumber

 typedef Vec_TNumber::iterator Vec_TNumber_it;   ///< iterator in Vec_TNumber
 typedef Vec_TNumber::const_iterator c_Vec_TNumber_it;
                                           ///< const iterator in Vec_TNumber

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 typedef double CNumber;                     ///< type of cost 
 typedef const CNumber c_CNumber;            ///< a read-only CNumber

 typedef std::vector< CNumber > Vec_CNumber;  ///< a vector of CNumber
 typedef const Vec_CNumber c_Vec_CNumber;     ///< a const vector of CNumber

 typedef Vec_CNumber::iterator Vec_CNumber_it;   ///< iterator in Vec_CNumber
 typedef Vec_CNumber::const_iterator c_Vec_CNumber_it;
                                           ///< const iterator in Vec_CNumber

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 typedef double TONumber; 
 /**< type of the objective function: has to hold sums of products of
    TNumber(s) by CNumber(s) */

 typedef const TONumber c_TONumber;             ///< a read-only TONumber

 typedef std::vector< TONumber > Vec_TONumber;  ///< a vector of TONumber
 typedef const Vec_TONumber c_Vec_TONumber;     ///< a const vector of TONumber

/** @} ---------------------------------------------------------------------*/
/*------------------------------- FRIENDS ----------------------------------*/
/*--------------------------------------------------------------------------*/

 friend OTSolution;  ///< make OTSolution friend

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of OTBlock, taking a pointer to the father (generic) Block
 /** Constructor of OTBlock. It accepts a pointer to the father Block, which
  * can be of any type, defaulting to nullptr so that this can also be used as
  * the void constructor. TODO: change below */

 explicit OTBlock( Block *father = nullptr ) 
  : Block( father ) , NAtomsP( 0 ) , NAtomsQ( 0 ) , NTransPlan( 0 ) , 
    MaxNAtomsP( 0 ) , MaxNAtomsQ( 0 ) ,  NStaticAtomsP( 0 ) , 
    NStaticAtomsQ( 0 ) , NStaticTransPlan( 0 ) , AR( 0 ) , f_cond_lower( 0 ) , 
    f_cond_upper( - Inf< double >() ) { }

/*--------------------------------------------------------------------------*/
 /// destructor of OTBlock: deletes the abstract representation, if any

 virtual ~OTBlock() { guts_of_destructor(); }

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// loads the OT instance from memory
 /** Loads the OT instance from memory. The parameters are what you expect:
  *
  * - n       number of atoms in probability distribution P
  * 
  * - m       number of atoms in probability distribution Q
  *
  * - wP      vector of proabibility weights of the probability distribution P
  *
  * - wQ      vector of probability weights of the probability distribution Q
  *
  * - cost    vector of unit costs of moving an unit i of the probability 
  *           distribution P to an unit j of the probability distribution Q.
  * 
  * - dynP    number (<= n and default 0) of dynamic atoms of P
  * 
  * - dynQ    number (<= m and default 0) of dynamic atoms of Q
  * 
  * - dynTP   number (<= n*m and default 0) of dynamic components of trans. plan
  * 
  * - mDynP   maximal number of dynamic atoms of P 
  * 
  * - mDynQ   maximal number of dynamic atoms of Q
  * 
  * TODO: add dynamic components post Monday discussion surely some atoms of 
  * Q should be put dynamic
  * 
  * */

 void load( Index n , Index m , c_Vec_TNumber& wP = {} , c_Vec_TNumber& wQ = {},
      c_Vec_CNumber & cost = {} , Index dynP = 0 , 
      Index dynQ = 0 , Index dynTP = 0 , Index mDynP = 0 , Index mDynQ = 0);

/*--------------------------------------------------------------------------*/
 /// load instance from txt file  
 /** Loads a OTBlock out of std::istream. The format is the
  * following, with each element being separated by whitespaces and possibly
  * comments:
  *
  * - n       number of atoms in probability distribution P
  * 
  * - m       number of atoms in probability distribution Q
  *
  * - wP      vector of proabibility weights of the probability distribution P
  *
  * - wQ      vector of probability weights of the probability distribution Q
  *
  * - cost    vector of unit costs of moving an unit i of the probability 
  *           distribution P to an unit j of the probability distribution Q.
  * 
  * - dynP    number (<= n and default 0) of dynamic atoms of P
  * 
  * - dynQ    number (<= m and default 0) of dynamic atoms of Q
  * 
  * - dynTP   number (<= n*m and default 0) of dynamic components of trans. plan
  * 
  * - mDynP   maximal number of dynamic atoms of P 
  * 
  * - mDynQ   maximal number of dynamic atoms of Q
  *
  * If the stream (after having extracted whitespaces and comments) does not
  * eof() here, then
  *
  * - for i = 1 to n: integrality of item i (true if integral, false if not)
  *
  * If integrality is not specified, true is assumed for all objects.
  *
  * Since there is only one supported input format, \p frmt is ignored.
  *
  * Like load( memory ), if there is any Solver attached to this 
  * OTBlock then a NBModification (the "nuclear option") is 
  * issued. */

 void load( std::istream & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// extends Block::deserialize( netCDF::NcGroup )
 /** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
  * a OTBlock. Besides what is managed by the serialize() method of the base
  * Block class, the group should contain the following:
  *
  *  - the dimension "NP" containing the number of atoms in the probability 
  *  distribution P;
  * 
  *  - the dimension "NQ" containing the number of atoms in the probability 
  *  distribution Q;
  * 
  *  - the variable "wP" of type double and indexed over the dimension "NP";
  *  the i-th component of "WP" contains the weight of the i-th atom of P.
  * 
  *  - the variable "wQ" of type double and indexed over the dimension "NQ";
  *  the j-th component of "wQ" contains the weight of the j-th atom of Q.
  * 
  *  - the variable "cost" of type double and indexed over the two dimensions
  *  "NP" and "NQ"; the (i*j + j)-th entry of the variable is assumed to contain 
  *  the cost of moving the i-th atom of P to the j-th atom of Q. 
  * 
  * All dimensions and variables are mandatory.  */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the OT problem
 /** Method that generates the abstract Variable of the OT problem. These are:    
  * 
  * TODO: */

 void generate_abstract_variables( Configuration *stvv = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the static constraint of the OT problem
 /** Method that generates the abstract constraint of the OT problem. These are:
  * TODO:
 */
 
 void generate_abstract_constraints( Configuration *stcc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the objective of the OT
 /** Method that generates the objective of the OT. 
  * 
  * *REMARK* On MCFBlock there is the possibility of alternating between sparse
  * and dense objective function. It could definitely be relevant for OTBlock as
  * the solution of the (unregularized) OT problem is sparse. However, it is
  * unclear if during an iterative algorithm that aims to compute the solution
  * of the OT problem the current transport plan is sparse or remains sparse.
  * 
  * TODO: Look into how to add optional regularization in the config?
  */

 void generate_objective( Configuration *objc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*------------ Methods for reading the data of the OTBlock -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of the OTBlock
 *  @{ */

 /// getting the current sense of the Objective, which is minimization

 [[nodiscard]] int get_objective_sense( void ) const override {
  return( Objective::eMin );
  }
  
/*--------------------------------------------------------------------------*/
 /// getting upper bounds on the value of the Objective
 /** An upper bound on the optimal value of the problem is computed as
  * \f[
  *  \sum_{ (i,j) \in N x M : c_{ i, j } > 0 } c_{ i, j } 
  * \f]
  * If it is finite (which it may not be), this is a valid upper bound.
  * */

 [[nodiscard]] double get_valid_upper_bound() {
  if( std::isnan( f_cond_upper ) )
   compute_conditional_bounds();

  return( f_cond_upper );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// getting a global valid lower bound on the value of the Objective
 /** A lower bound on the optimal value of the problem is 0. A tighter lower 
  * bound can be computed as
  * \f[
  *  \sum_{ j \in M } [ min_{ i \in N } c_{ i, j } ].
  * \f]
  * That is, for every atom j of Q, we send all its mass to the atom i of P 
  * which has the lowest cost. It might not be associated to an admissible 
  * transport plan (as the atom i might not be able to hold all the mass of j)
  * but it yields a lower bound on the optimal value nonetheless. */

 [[nodiscard]] double get_valid_lower_bound()
  {
    return( f_cond_lower );
  }

/*--------------------------------------------------------------------------*/
 /// get the number of atoms in P

 [[nodiscard]] Index get_NAtomsP( void ) const { return( NAtomsP ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the number of atoms in Q
 [[nodiscard]] Index get_NAtomsQ( void ) const { return( NAtomsQ ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the number of static arcs

 [[nodiscard]] Index get_NTransPlan( void ) const 
{
 return( NTransPlan );
}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the maximum number of atoms in P

 [[nodiscard]] Index get_MaxNAtomsP( void ) const { return( MaxNAtomsP ); }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the maximum number of atoms in Q

 [[nodiscard]] Index get_MaxNAtomsQ( void ) const { return( MaxNAtomsQ ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the maximum number of components in the transport plan

 [[nodiscard]] Index get_MaxNTransPlan( void ) const { 
  return( MaxNAtomsP*MaxNAtomsQ ); 
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the number of static atoms in P

 [[nodiscard]] Index get_NStaticAtomsP( void ) const {
  return( NStaticAtomsP );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the number of static atoms in Q

 [[nodiscard]] Index get_NStaticAtomsQ( void ) const {
  return( NStaticAtomsQ );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the number of static arcs

 [[nodiscard]] Index get_NStaticTransPlan( void ) const {
  return( NStaticTransPlan );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if there are static atoms in P

 [[nodiscard]] bool HasStaticP( void ) const {
  return( get_NStaticAtomsP() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if there are static atoms in Q

 [[nodiscard]] bool HasStaticQ( void ) const {
  return( get_NStaticAtomsQ() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if there are static transport plan components

 [[nodiscard]] bool HasStaticTransPlan( void ) const { 
  return( get_NStaticTransPlan() ); 
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if there are dynamic atoms in P

 [[nodiscard]] bool HasDynamicP( void ) const {
  return( get_NAtomsP() > get_NStaticAtomsP() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if there are dynamic atoms in Q

 [[nodiscard]] bool HasDynamicQ( void ) const {
  return( get_NAtomsQ() > get_NStaticAtomsQ() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if there are dynamic transport plan components

 [[nodiscard]] bool HasDynamicTransPlan( void ) const {
  return( get_NTransPlan() > get_NStaticTransPlan() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if there may ever be dynamic atoms in P

 [[nodiscard]] bool MayHaveDynP( void ) const {
  return( get_MaxNAtomsP() > get_NStaticAtomsP() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if there may ever be dynamic atoms in Q

 [[nodiscard]] bool MayHaveDynQ( void ) const {
  return( get_MaxNAtomsQ() > get_NStaticAtomsQ() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if there may ever be dynamic transport plan components

 [[nodiscard]] bool MayHaveDynTransPlan( void ) const {
  return( get_MaxNTransPlan() > get_NStaticTransPlan() );
  }
  
/*--------------------------------------------------------------------------*/
 /// given a pointer to a transport plan Variable, return component index
 /** Given a pointer to a trans. plan Variable (formally a Variable *, but
  * immediately static_cast-ed to a ColVariable * right inside), returns the
  * index of the corresponding component. Throws exception if the pointer is 
  * not to a [Col]Variable of the OTBlock. */

 [[nodiscard]] Index p2i_x( const Variable * var ) const {
  auto i = p2i_x_s( var );
  if( ( i >= 0 ) && ( i < int( get_NStaticTransPlan() ) ) )
   return( i );

  i = get_NStaticTransPlan();
  for( auto dxi = dx.begin() ; dxi != dx.end() ; ++i , ++dxi )
   if( &(*dxi) == static_cast< const ColVariable * >( var ) )
    return( i );

  throw( std::invalid_argument( "invalid transport plan pointer" ) );
  return( 0 );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// given a transport plan component, return the pointer to trans plan variable
 /** Given the index of an arc, returns the pointer to the corresponding 
  * transport plan variable (a ColVariable *). 
  * This ASSUMES THE Variable ARE CONSTRUCTED IN THE FIRST PLACE, SEGFAULTS ARE 
  * BOUND TO HAPPEN OTHERWISE. */

 [[nodiscard]] ColVariable * i2p_x( Index i ) const {
  if( i >= get_NTransPlan() )
   throw( std::invalid_argument( "invalid transport plan name" ) );

  if( i < get_NStaticTransPlan() )
   return( const_cast< ColVariable * >( & x[ i ] ) );
  else
   if( i - get_NStaticTransPlan() < get_NTransPlan() - i )
    return( const_cast< ColVariable * >(
		   &( *std::next( dx.begin() , i - get_NStaticTransPlan() ) ) ) );
   else
    return( const_cast< ColVariable * >(
		           &( *std::prev( dx.end() , get_NTransPlan() - i ) ) ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if the component is closed; deleted component are not closed

 [[nodiscard]] bool is_closed( Index component ) const {
  return( ( ! is_deleted( component ) ) && i2p_x( component )->is_fixed() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if the component is deleted

 [[nodiscard]] bool is_deleted( Index component ) const {
  return( ( ! vecC.empty() ) && std::isnan( vecC[ component ] ) );
  }
  
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the vector of arc costs
 /** Returns a const reference to the vector of arc costs of size
  * get_MaxNTransPlan(). Note that the cost of a deleted arc is NaN. */

 [[nodiscard]] c_Vec_CNumber & get_vecC( void ) const { return( vecC ); }

 /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// get the cost of arc i (0 <= i < get_NTransPlan()), NaN if deleted

 [[nodiscard]] CNumber get_vecC( c_Index i ) const { return( vecC[ i ] ); }

/** @} ---------------------------------------------------------------------*/
/*--------------------- Methods for checking the Block ---------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking the Block
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// returns true if the current solution is approximately feasible
 /** Returns true if the solution encoded in the current value of the transport
  * plan (x) Variable of the OTBlock is approximately feasible. This clearly
  * requires the Variable of the OTBlock to have been defined, i.e., that
  * generate_abstract_variables() has been called prior to this method.
  *
  * The parameter for deciding what "approximately feasible" exactly means is
  * a single TNumber value, representing the *relative* tolerance for
  * satisfaction of both flow conservation constraint and flow upper/lower
  * bounds. This value is to be found as:
  *
  * - if fsbc is not nullptr and it is a SimpleConfiguration< FNumber >, then
  *   it is fsbc->f_value;
  *
  * - otherwise, if f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_is_feasible_Configuration is not nullptr and it
  *   is a SimpleConfiguration< FNumber >, then it is
  *   f_BlockConfig->f_is_feasible_Configuration->f_value;
  *
  * - otherwise, it is 0. */
 
 bool is_feasible( bool useabstract = false , Configuration *fsbc = nullptr )
  override;

/** @} ---------------------------------------------------------------------*/
/*------------------------- Methods for R3 Blocks --------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for R3 Blocks
 *  @{ */

 /// gets an R3 Block of OTBlock currently only the copy one
 /** Gets an R3 Block of the OTBlock. The list of currently supported R3
  * Block is:
  *
  * - r3bc == nullptr: the copy (an OTBlock identical to this)
  */

 Block * get_R3_Block( Configuration *r3bc = nullptr ,
		       Block * base = nullptr , Block * father = nullptr )
  override;

/*--------------------------------------------------------------------------*/
 /// maps back the solution from a copy OTBlock to the current one
 /** Maps back the solution from a copy OTBlock to the current one. The
  * parameter r3bc is useless (has to be nullptr). The parameter solc decides
  * which part of the solution is mapped:
  *
  * - if solc != nullptr and it is a SimpleConfiguration< int >, then it
  *   depends on solc->f_value:
  *
  *   = 1 means "only map the primal solution"
  *
  *   = 2 means "only map the dual solution"
  *
  *   = everything else (e.g., 0) means "map everything";
  *
  * - if solc == nullptr, f_BlockConfig != nullptr,
  *   f_BlockConfig->f_solution_Configuration != nullptr and it
  *   is a SimpleConfiguration< int >, then it depends on its f_value as in
  *   the previous case;
  *
  * - otherwise, everything (both the primal and the dual solution) is
  *   mapped.
  *
  * The same format applies verbatim to the case of primal or dual unbounded
  * rays (negative-cost unbounded cycles and cuts, respectively), although
  * one would expect only one of these to be found (but both may
  * theoretically do).
  *
  * Note that R3B may not contain some or all of the required solution, if
  * the corresponding Variable/Constraint have not been constructed yet:
  * this throws an exception. */ 

 void map_back_solution( Block *R3B , Configuration *r3bc = nullptr ,
			 Configuration *solc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// maps the solution of the current OTBlock to a copy OTBlock
 /** Maps the solution of the current OTBlock to a copy OTBlock. The
  * parameter r3bc is useless (has to be nullptr). The parameter solc decides
  * which part of the solution is mapped:
  *
  * - if solc != nullptr and it is a SimpleConfiguration< int >, then it
  *   depends on solc->f_value:
  *
  *   = 1 means "only map the primal solution"
  *
  *   = 2 means "only map the dual solution"
  *
  *   = everything else (e.g., 0) means "map everything";
  *
  * - if solc == nullptr, f_BlockConfig != nullptr,
  *   f_BlockConfig->f_is_solution_Configuration != nullptr and it
  *   is a SimpleConfiguration< int >, then it depends on its f_value as in
  *   the previous case;
  *
  * - otherwise, everything (both the primal and the dual solution) is
  *   mapped.
  *
  * The same format applies verbatim to the case of primal or dual unbounded
  * rays (negative-cost unbounded cycles and cuts, respectively), although
  * one would expect only one of these to be found (but both may
  * theoretically do).
  *
  * Note that the current OTBlock may not contain some or all of the
  * required solution, if the corresponding Variable/Constraint have not
  * been constructed yet: this throws an exception. */ 

 void map_forward_solution( Block *R3B , Configuration *r3bc = nullptr ,
			    Configuration *solc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /** No specific Configuration is required, hence expected, for OTBlock.
  *
  * IMPORTANT NOTE: map_forward_Modification() only maps "physical"
  * Modification. The point is that if any part of the "abstract
  * representation" of OTBlock is changed, the corresponding "abstract"
  * Modification is intercepted in add_Modification() and a "physical"
  * Modification is also issued. Hence, for any change in OTBlock there
  * will always be both Modification "in flight", and therefore there is
  * no need (and good reasons not) to map both.
  *
  * In particular, the method handles the following Modification:
  *
  * - GroupModification
  *
  * - OTBlockRngdMod
  *
  * - OTBlockSbstMod
  *
  * - NBModification
  *
  * Any other Modification is ignored (and false is returned).
  *
  *     IMPORTANT NOTE: OTBlockRngdMod ALLOW TO ADD/DELETE ARCS IN THE
  *     PROBLEM, WHICH ALSO CHANGES THE "NAMES" OF EXISTING ARCS. OTBlock
  *     IMPLEMENTS map_forward_Modification() IN A WAY THAT IS ONLY
  *     GUARANTEED TO BE CORRECT IF:
  *
  *     = EITHER THE SET OF ARCS IS NEVER CHANGED;
  *
  *     = OR THE Modification ARE MAPPED IMMEDIATELY AFTER THEY ARE ISSUED.
  *
  * This is because otherwise OTBlock should have to understand whether the
  * set of arc "names" in the Modification is still correct and do something
  * in case it is not, which is too complex to do at the moment.
  *
  * Note that for GroupModification, true is returned only if all the
  * inner Modification of the GroupModification return true.
  *
  * Note that if the issueAMod param is eModBlck, then it is "downgraded" to
  * eNoBlck: the method directly does "physical" changes, hence there is no
  * reason for it to issue "abstract" Modification with concerns_Block() ==
  * true. */

 bool map_forward_Modification( Block *R3B , c_p_Mod mod ,
				Configuration *r3bc = nullptr ,
				ModParam issuePMod = eNoBlck ,
				ModParam issueAMod = eModBlck ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /** No specific Configuration is required, hence expected, for OTBlock.
  *
  * The current implementation of map_back_Modification() actually uses
  * map_forward_Modification() in reverse, so see the comments to the latter
  * method. */

 bool map_back_Modification( Block *R3B , c_p_Mod mod ,
			     Configuration *r3bc = nullptr ,
			     ModParam issuePMod = eNoBlck ,
			     ModParam issueAMod = eModBlck ) override;

/** @} ---------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 *  @{ */

 /// returns a OTSolution representing the current solution of this OTBlock
 /** Returns a OTSolution representing the current solution status of this
  * OTBlock. What kind of solution is saved depends on the integer value ws,
  * obtained as follows:
  *
  * - if solc != nullptr and it is a SimpleConfiguration< int >, then
  *   ws == solc->f_value:
  *
  * - if solc == nullptr, f_BlockConfig != nullptr,
  *   f_BlockConfig->f_solution_Configuration != nullptr and it
  *   is a SimpleConfiguration< int >, ws is its f_value
  *
  * - otherwise ws is 0.
  *
  * The encoding of ws is:
  *
  *   = 1 means "only save the primal solution"
  *
  *   = 2 means "only save the dual solution"
  *
  *   = everything else (e.g., 0) means "save everything";
  *
  * The same format applies verbatim to the case of primal or dual unbounded
  * rays (negative-cost unbounded cycles and cuts, respectively), although
  * one would expect only one of these to be found (but both may
  * theoretically do).
  *
  * Note that OTBlock may not contain some or all of the required solution,
  * if the corresponding Variable/Constraint have not been constructed yet:
  * this throws an exception, unless emptys = true, in which case the
  * OTSolution object is only prepped for getting a solution, but it is not
  * really getting one now.
  *
  * Note that, although the method clearly returns a OTSolution, formally
  * the return type is Solution *. This is because it is not possible to
  * forward declare OTSolution as a derived class from Solution, nor to
  * define OTSolution before OTBlock because the former uses some type
  * information declared in the latter. */ 

 Solution * get_Solution( Configuration *solc = nullptr ,
			  bool emptys = true ) override;

 /*--------------------------------------------------------------------------*/
 /// returns the objective value of the current solution

 TONumber get_objective_value( void ) {
  if( ! ( AR & HasObj ) )  // the objective is not there
   return( Inf< RealObjective::OFValue >() );
  c.compute();
  return( c.value() );
  }

/*--------------------------------------------------------------------------*/
 /// gets a contiguous interval of the flow solution
 /** Method to get the flow solution; upon return, the current value of the
  * flow solution for the i-th arc in \p rng is written in *( FSol + i ).
  * Note that if the right extreme of the range is >= get_NTransPlan() it is
  * ignored. */

 void get_x( Vec_TNumber_it FSol , Range rng = Range( 0 , Inf< Index >() ) )
  const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// gets the flow solution for an arbitrary subset of arcs
 /** Method to get the flow solution; upon return, the current value of the
  * flow solution for arc nms[ i ] for all 0 <= i <  nms.size() is written
  * in *( FSol + i ). Note that
  *
  *     nms IS ASSUMED TO BE ORDERED BY INCREASING Index */

 void get_x( Vec_TNumber_it FSol , c_Subset & nms ) const;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// gets the transport plan component value of a cost arc

 TNumber get_x( Index arc ) const {
  if( arc >= get_NTransPlan() )
   throw( std::invalid_argument( "invalid cost arc name" ) );

  if( arc < get_NStaticTransPlan() )
   return( x[ arc ].get_value() );
  else
   return( std::next( dx.begin() , arc - get_NStaticTransPlan() )->get_value() );
  }

/*--------------------------------------------------------------------------*/
 /// gets a contiguous interval of the potential solution
 /** Method to get the potential solution; upon return, the current value of
  * the potential solution for the i-th node in \p rng is written into
  * *( PSol + i ). Note that if the right extreme of the range is >=
  * get_NNodes() it is ignored. Note that "node names" here go from 0 to
  * get_NNodes() - 1, despite the fact that get_STP() and get_EN() report node
  * "names" between 1 and get_NNodes(). */

 void get_pi( Vec_CNumber_it PSol , Range rng = Range( 0 , Inf< Index >() ) )
  const;

/*--------------------------------------------------------------------------*/
 /// sets a contiguous interval of the transport plan components
 /** Method to set the transport plan components; the values found in the 
  * c_Vec_TNumber starting from fstrt are copied into the value of the transport
  * plan variable x[ i ] for i in rng, in the same order. */

 void set_x( c_Vec_TNumber_it fstrt ,
	     Range rng = Range( 0 , Inf< Index >() ) );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets a generic subset of the transport plan component
 /** Method to set the transport plan component; the values found in the 
  * c_Vec_TNumber starting from fstrt are copied into the value of the flow 
  * variable x[ i ] for all i in sbst (that must be ordered in increasing 
  * sense), in the same order. */

 void set_x( c_Vec_TNumber_it fstrt , c_Subset sbst );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the transport plan solution of the given arc

 void set_x( Index arc , TNumber FSol ) {
  if( arc >= get_NTransPlan() )
   throw( std::invalid_argument( "invalid arc name" ) );

  if( arc < get_NStaticTransPlan() )
   x[ arc ].set_value( FSol );
  else
   std::next( dx.begin() , arc - get_NStaticTransPlan() )->set_value( FSol );
  }

/** @} ---------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Modification
 *  @{ */

 /// returns true if there is any Solver "listening to this OTBlock"
 /** Returns true if there is any Solver "listening to this OTBlock", or if
  * the OTBlock has to "listen" anyway because the "abstract" representation
  * is constructed, and therefore "abstract" Modification have to be generated
  * anyway to keep the two representations in sync.
  *
  * No, this should not be needed. In fact, if the "abstract" representation
  * is modified with the default eModBlck value of issueMod, it is issued
  * irrespectively to the value of anyone_there(); see Observer::issue_mod().
  * If the value of issueMod is anything else the  "abstract" representation
  * has been modified already and there is no point in issuing the
  * Modification.
  * Note that that Observer::issue_mod() does not check if the "abstract"
  * representation has been constructed, but this is clearly not
  * necessary, as the Modification we are speaking of are issued while
  * changing the "abstract" representation, if that has not been
  * constructed then it cannot issue Modification

 bool anyone_there( void ) const override {
  return( AR ? true : Block::anyone_there() );
  }
 */
/*--------------------------------------------------------------------------*/
 /// adding a new Modification to the OTBlock
 /** Method for handling Modification.
  *
  * The version of OTBlock has to intercept any "abstract Modification" that
  * modifies the "abstract representation" of the OTBlock, and "translate"
  * them into both changes of the actual data structures and corresponding
  * "physical Modification". These Modification are those for which
  * Modification::concerns_Block() is true. Note, however, that before sending
  * the Modification to the Solver and/or the father Block, the
  * concerns_Block() value is set to false. This is because once it is passed
  * through this method, the "abstract Modification" has "already done its
  * duty" of providing the information to the OTBlock, and this must not be
  * repeated. In particular, this would be an issue if the Modification would
  * be [map_forward or map_back]-ed, because inside of this method a "physical
  * Modification" doing the same job is surely issued. That Modification would
  * also be [map_forward or map_back]-ed, together with the original "abstract
  * Modification" that would pass again through this method (in the other
  * OTBlock), which would mean that the "physical Modification" would be
  * issued twice.
  *
  * The following "abstract Modification" are handled:
  *
  * - GroupModification, that are simply unpacked into the individual
  *   sub-[Group]Modification and dealt with individually;
  *
  * - C05FunctionModRngd and C05FunctionModSbst changing coefficients coming
  *   from the (LinearFunction into the FRow)Objective, but *not* from the
  *   (LinearFunction into the FRow)Constraint;
  *
  * - RowConstraintMod changing the RHS of the bound constraints and both
  *   sides at once of the flow conservation ones, but not any other
  *   combination; and note that the RHS of the bound constraints may not
  *   be changeable at all if they have not been constructed, in which
  *   case there cannot be any Modification to handle here;
  *
  * - VariableMod fixing and un-fixing a flow ColVariable; however, note
  *   that *fixing is only permitted if the value() of the ColVariable is
  *   zero*, because that corresponds to closing the arc, exception being
  *   thrown otherwise.
  *
  * Any other Modification reaching the OTBlock will lead to exception
  * being thrown.
  *
  * Note: any "physical" Modification resulting from processing an "abstract"
  *       one will be sent to the same channel (chnl). */

 void add_Modification( sp_Mod mod , ChnlName chnl = 0 ) override;

/** @} ---------------------------------------------------------------------*/
/*--------------- METHODS FOR PRINTING & SAVING THE OTBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the OTBlock
 *  @{ */

 /// print the OTBlock on an ostream with the given verbosity
 /** Protected method to print information about the OTBlock; with the
  * "complete" level ('C') it outputs the OTBlock in DIMACS format. */

 void print( std::ostream & output , char vlvl = 0 ) const override;

/*--------------------------------------------------------------------------*/
 /// extends Block::serialize( netCDF::NcGroup )
 /** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
  * OTBlock. See OTBlock::deserialize( netCDF::NcGroup ) for details of the
  * format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*------------- METHODS FOR ADDING / REMOVING / CHANGING DATA --------------*/
/*--------------------------------------------------------------------------*/
/** @name Changing the data of the OT instance
 *
 * All the methods in this section have two parameters issueMod and issueAMod
 * which control if and how the, respectively, "physical Modification" and
 * "abstract Modification" corresponding to the change have to be issued, and
 * where (to which channel). The format of the parameters is that of
 * Observer::make_par(), except that the value eModBlck is ignored and
 * treated it as if it were eNoBlck [see Observer::issue_pmod()]. This is
 * because it makes no sense to issue an "abstract" Modification with
 * concerns_Block() == true, since the changes in the OTBlock have surely
 * been done already, and this is just not possible for a "physical"
 * Modification.
 *
 * IMPORTANT NOTE: the current implementation of all these methods issues (at
 * most) *two separate* Modification, a "physical" and an "abstract" one. The
 * latter may be a GroupModification bunching together related abstract
 * Modification, but the two Modification are nonetheless separate. A
 * different approach could be to issue a single GroupModification with inside
 * both the "physical" and the "abstract" one (the latter possibly itself a
 * GroupModification). This may allow a more efficient handling of
 * Modification by ensuring that the two are always received together, but at
 * the cost of a more intricate code that is best avoided for now.
 *
 * Note: the methods accept the eDryRun value for the issueAMod parameter for
 * the "abstract" representation. This allows to re-use them within OTBlock
 * itself when reacting to abstract Modification, where the  "abstract"
 * representation has been changed already.
 *  @{ */

 /// change the costs of a contiguous interval of arcs
 /** Method to change the costs of a subset of arcs with "contiguous names".
  * That is, *( NCost + i - strt ) becomes the new cost of the i-th arc in
  * \p rng. Note that if the right extreme of the range is >= get_NTransPlan() it 
  * is ignored.
  *
  * Note that if \p rng contains some closed arc, its cost is also changed.
  * While this has no immediate impact on the problem solved, if the arc is
  * re-opened then the cost set with this method when the arc was closed is
  * in effect.
  *
  * If more than one Modification is actually issued and issueAMod specifies
  * an open channel, then the channel is nested so that the three Modification
  * are grouped into a single GroupModification. Similarly, if instead
  * issueAMod specifies the default channel, then a new channel is opened to
  * group the multiple Modification and immediately closed when the last one
  * is issued. If, instead, the Objective is a "dense" LinearFunction, then
  * at most one LinearFunctionMod for modifying the coefficients is issued.
  * Of course this only applies if issueAMod specifies that abstract
  * Modification have to be issued *and* the abstract Objective has been
  * constructed.
  *
  * Also, if issueMod says so then a "physical" OTBlockRngdMod is issued. */

 void chg_costs( c_Vec_CNumber_it NCost , Range rng = INFRange ,
		 ModParam issueMod = eNoBlck , ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// change the costs of an arbitrary subset of arcs
 /** Method to change the costs of an arbitrary subset of arc. That is,
  * *( NCost + i ) becomes the new cost of arc nms[ i ] for all 0 <= i <
  * NCost.size(), (which means that nms.size() == NCost.size()). The
  * parameter ordered tells if the nms vector is ordered for increasing
  * index of the arc. As the && tells, nms is "consumed" by the method,
  * typically being shipped to an appropriate OTBlockSbstMod object.
  *
  * See chg_costs( range ) for Modification issued (except that, of course,
  * the "physical" one is a OTBlockSbstMod), and about changes in costs
  * of closed arcs. */

 void chg_costs( c_Vec_CNumber_it NCost ,
		 Subset && nms , bool ordered = false ,
		 ModParam issueMod = eNoBlck , ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// changes the cost of the given arc
 /** Changes the cost of the given arc.
  *
  * Note that this can issue only one Modification of each type; the
  * "physical" one is a OTBlockRngdMod with rng = [ arc ). */

 void chg_cost( CNumber NCost , Index arc ,
		ModParam issueMod = eNoBlck , ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// closes a contiguous interval of arcs
 /** Method to close a subset of arcs with all "contiguous names" given in
  * \rng; note that if the right extreme of the range is >= get_NTransPlan() it is
  * ignored. The flow on the arcs is fixed to 0 but the arcs are not removed
  * from the problem, and their capacity and cost are not changed, so that
  * they can be easily re-opened later. When the problem is created, all arcs
  * are open. Closing an already closed arc does nothing.
  *
  * Note that closing multiple arcs can issue as many Modification as there
  * are arcs in the range, in particular VariableMod. If more than one
  * Modification is actually issued and issueAMod specifies an open channel,
  * then the channel is nested so that all the Modification are grouped into
  * a single GroupModification. Similarly, if instead issueAMod specifies the
  * default channel, then a new channel is opened to group the multiple
  * Modification and immediately closed when the last one is issued. Of course
  * this only applies if issueAMod specifies that abstract Modification have
  * to be issued.
  *
  * Also, if issueMod says so then a "physical" OTBlockRngdMod is issued. */

 void close_arcs( Range rng = INFRange ,
		  ModParam issueMod = eNoBlck ,
		  ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// closes an arbitrary subset of arcs
 /** Method to close an arbitrary subset of arc, i.e., all those whose names
  * are found in the array nms. The flow on the arcs is fixed to 0 but the
  * arcs are not removed from the problem, and their capacity and cost are
  * not changed, so that they can be easily re-opened later. When the problem
  * is created, all arcs are open. Closing an already closed arc does
  * nothing.
  *
  * The parameter ordered tells if the nms vector is ordered for increasing 
  * index of the arc. As the && tells, nms is "consumed" by the method,
  * typically being shipped to an appropriate OTBlockSbstMod object.
  *
  * See close_arcs( range ) for Modification issued (except that, of course,
  * the "physical" one is a OTBlockSbstMod). */

 void close_arcs( Subset && nms , bool ordered = false ,
		  ModParam issueMod = eNoBlck ,
		  ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// closes the given arc
 /** Method to "close" the given arc: the flow on arc is fixed to 0. The arc
  * is not removed from the problem, and its capacity and cost are not
  * changed, so that it can be easily re-opened later. When the problem is
  * created, all arcs are open. Closing an already closed arc does nothing.
  *
  * Note that this can issue only one Modification; the "physical" one is a
  * OTBlockRngdMod with rng = [ arc ). */

 void close_arc( Index arc , ModParam issueMod = eNoBlck ,
		             ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
/// re-opens a contiguous interval of arcs
 /** Method to "open" a subset of closed arcs with all "contiguous names"
  * given in \rng; note that if the right extreme of the range is >=
  * get_NTransPlan() it is ignored. Opening an already open arc (which is what all
  * arcs are when the problem is created) does nothing.
  *
  * Note that opening multiple arcs can issue as many Modification as there
  * are arcs in the range, in particular VariableMod. If more than one
  * Modification is actually issued and issueAMod specifies an open channel,
  * then the channel is nested so that all the Modification are grouped into
  * a single GroupModification. Similarly, if instead issueAMod specifies the
  * default channel, then a new channel is opened to group the multiple
  * Modification and immediately closed when the last one is issued. Of course
  * this only applies if issueAMod specifies that abstract Modification.
  *
  * Also, if issueMod says so then a "physical" OTBlockRngdMod is issued. */

 void open_arcs( Range rng = INFRange ,
		 ModParam issueMod = eNoBlck ,
		 ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// re-opens an arbitrary subset of arcs
 /** Method to "open" an arbitrary subset of closed arc, i.e., all those
  * whose names are found in the array nms. Opening an already open arc
  * (which is what all arcs are when the problem is created) does nothing.
  *
  * The parameter ordered tells if the nms vector is ordered for increasing 
  * index of the arc. As the && tells, nms is "consumed" by the method,
  * typically being shipped to an appropriate OTBlockSbstMod object.
  *
  * Note that closing multiple arcs can issue as many Modification as there
  * are arcs in the range, in particular VariableMod. If more than one
  * Modification is actually issued and issueAMod specifies an open channel,
  * then the channel is nested so that all the Modification are grouped into
  * a single GroupModification. Similarly, if instead issueAMod specifies the
  * default channel, then a new channel is opened to group the multiple
  * Modification and immediately closed when the last one is issued.
  * Of course this only applies if issueAMod specifies that abstract
  * Modification have to be issued *and* the abstract Constraint have been
  * constructed.
  *
  * Also, if issueMod says so then a "physical" OTBlockRngdMod is issued. */

 void open_arcs( Subset && nms , bool ordered = false ,
		 ModParam issueMod = eNoBlck ,
		 ModParam issueAMod = eNoBlck );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// re-opens the given arc
 /** Method to "open" the given closed arc, i.e., allow the flow on arc to
  * vary. Opening an already open arc (which is what all arcs are when the
  * problem is created) does nothing.
  *
  * Note that this can issue only one Modification; the "physical" one is a
  * OTBlockRngdMod with rng = [ arc ). */

 void open_arc( Index arc , ModParam issueMod = eNoBlck ,
		            ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// add a new arc
 /** Method to add a new arc, providing its starting and ending nodes, cost
  * and capacity.
  *
  * The method returns the "name" that the new arc has received, which is the
  * index by which the arc has to be addressed in all methods (such as
  * chg_[cost/ucap](), [open/close]_arc(), get_[x/rc]()). The chosen name
  * depends on whether or not there are "deleted" arcs (see remove_arc())
  * with "name" < get_NTransPlan() - 1. If there is any such arc, then the "name"
  * of the new arc will be the smallest index among these; in this case, the
  * value returned by get_NTransPlan() does *not* change. Otherwise, if get_NTransPlan()
  * < get_MaxNTransPlan() then the new arc gets "name" get_NTransPlan(), which is
  * returned by the method, and the value returned by get_NTransPlan() increases by
  * one. Otherwise the arc is not actually added, and the method returns
  * Inf< FNumber >().
  *
  * Successfully adding a new arc causes the issuing of several Modification,
  * unless the issueMod and issueAMod parameters prevent this to happen:
  *
  * - a "physical" OTBlockRngdMod with type eAddArc;
  *
  * - an "abstract" GroupModification containing up to:
  *
  *   = two C05FunctionModVars corresponding to having added the new Variable
  *     to the two flow conservation constraints of its starting and ending
  *     node;
  *
  *   = one OneVarConstraintMod with type RowConstraintMod::eChgRHS for
  *      modifying the flow bound;
  *
  *   = if the "name" of the arc is == get_NTransPlan() (before the call):
  *
  *     * a BlockModAdd< ColVariable > corresponding to the addition of a
  *       new dynamic Variable (the flow Variable of the arc);
  *
  *     * possibly, a BlockModAdd< LB0Constraint > corresponding to the
  *       addition of a new dynamic Constraint (the bound Constraint of the
  *       arc, if it is defined);
  *
  *     * one more C05FunctionModVars corresponding to having added the new
  *       Variable to the objective.
  *
  *   = if, instead, the "name" of the arc is < get_NTransPlan() (before the call):
  *
  *     * one C05FunctionModLin for modifying the cost coefficients;
  *
  *     * one VariableMod making the flow variable "free";
  *       
  *  Of course, all the "abstract" Modification are only issued if the
  *  corresponding part of the "abstract" representation is constructed. */
 
 Index add_arc( Index sn , Index en , CNumber cst = 0 ,
		TNumber cap = Inf< TNumber >() ,
		ModParam issueMod = eNoBlck , ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// removes an existing arc
 /** Method to remove the arc which given name. It must be
  * get_NStaticTransPlan() <= arc < get_NTransPlan(), otherwise exception is thrown.
  *
  * The operation is performed differently in the case where arc ==
  * get_NTransPlan() - 1, i.e., the very last arc is eliminated, or
  * arc < get_NTransPlan() - 1.
  *
  * In the latter case, the elimination of the arc is "virtual", in the
  * sense that the flow variable is kept, together with all corresponding
  * parts of the "abstract" representation (if constructed). Only, the
  * value of the flow Variable is changed to 0 and the Variable is fixed, as
  * when the arc is closed. Furthermore, its starting and ending nodes (as
  * returned by get_STP() and get_EN()) are set to Inf< Index >(). This means
  * that the value returned by get_NTransPlan() does *not* change.
  *
  * In the former case, the elimination of the arc is "physical": not only
  * of that arc, but also of and all the "deleted" arcs with smaller name
  * up until the first non-deleted arc (or get_NStaticTransPlan()). The value
  * of get_NTransPlan() changes accordingly, and all the corresponding parts of
  * the "abstract" representation (Variable, bound constraints, coefficients
  * in the objective function and the constraints) are removed.
  *
  * Removing an existing arc causes the issuing of several Modification,
  * unless the issueMod and issueAMod parameters prevent this to happen:
  *
  * - a "physical" OTBlockRngdMod with type eRmvArc;
  *
  * - an "abstract" GroupModification containing up to:
  *
  *   = for each of the removed arcs, two C05FunctionModVars corresponding
  *     to having removed the existing Variable from the two flow
  *     conservation constraints of its starting and ending node;
  *
  *   = if the elimination is "virtual", one VariableMod corresponding to
  *     fixing the flow variable;
  *
  *   = if the elimination is "physical":
  *
  *     * a BlockModRmv< ColVariable > corresponding to the removal of the
  *       existing dynamic Variable (the flow Variable of the arcs);
  *
  *     * possibly, a BlockModRmv< LB0Constraint > corresponding to the
  *       removal of the existing dynamic Constraint (the bound Constraint of
  *       the arcs, if they are defined);
  *
  *     * for each of the removed arcs, one more C05FunctionModVars
  *       corresponding to having removed the existing Variable from the
  *       objective.
  *
  *  Of course, all the "abstract" Modification are only issued if the
  *  corresponding part of the "abstract" representation is constructed. */

 void remove_arc( Index arc , ModParam issueMod = eNoBlck ,
		              ModParam issueAMod = eNoBlck );

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED FIELDS  ----------------------------*/
/*--------------------------------------------------------------------------*/

 Index NAtomsP;                  ///< the current number of atoms of P
 Index NAtomsQ;                  ///< the current number of atoms of Q
 Index NTransPlan;               ///< current number components of trans. plan.
 Index MaxNAtomsP;               ///< the maximum number of atoms of P
 Index MaxNAtomsQ;               ///< the maximum number of atoms of Q
 Index NStaticAtomsP;            ///< the number of static atoms of P
 Index NStaticAtomsQ;            ///< the number of static atoms of Q
 Index NStaticTransPlan;         ///< number static components of trans. plan.

 Vec_CNumber vecC;               ///< vector of arc costs
 Vec_TNumber U;                  ///< vector holding the transport plan 

 unsigned char AR;               ///< bit-wise coded: what abstract is there

 static constexpr unsigned char HasVar = 1;
 ///< first bit of AR == 1 if the Variable have been constructed
 static constexpr unsigned char HasObj = 2;
 ///< second bit of AR == 1 if the Objective has been constructed
 static constexpr unsigned char HasFlw = 4;
 ///< third bit of AR == 1 if the Flow Conservation have been constructed
 static constexpr unsigned char HasBnd = 8;
 ///< fourth bit of AR == 1 if the Bound have been constructed

 double f_cond_lower;            ///< conditional lower bound, can be -INF
 double f_cond_upper;            ///< conditional upper bound, can be +INF
 
 std::vector< ColVariable > staticTP;          ///< the static TP variables
 std::vector< FRowConstraint> massP;    ///< static TP mass constraints on P
 std::vector< FRowConstraint> massP;    ///< static TP mass constraints on Q
 std::vector< LB0Constraint > UB;       ///< the static bound constraints
 
 std::list< ColVariable > dynamicTP;    ///< the dynamic TP variables
 std::list< FRowConstraint > dmassP;    ///< dynamic mass constrs. on P
 std::list< FRowConstraint > dmassQ;    ///< dynamic mass constrs. on P

 std::list< LB0Constraint > dUB;        ///< the dynamic bound constraints

 FRealObjective objective;               ///< the (linear) objective function

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/// register OTBlock methods into the method factories
/** Although in general private methods should not be commented, this one is
 * because it does the registration of the following OTBlock methods:
 *
 * - chg_costs() (both range and subset version)
 *
 * - chg_ucaps() (both range and subset version)
 *
 * - chg_dfcts() (both range and subset version)
 *
 * - close_arcs() (both range and subset version)
 *
 * - open_arcs() (both range and subset version)
 *
 * into the corresponding method factories. */

 static void static_initialization( void )
 {
  /*!!
 * Not all C++ compilers enjoy the template wizardry behind the three-args
 * version of register_method<> with the compact MS_*_*::args(), so we just
 * use the slightly less compact one with the explicit argument and be done
 * with it. !!*/
  // register_method< OTBlock >( "OTBlock::chg_costs", &OTBlock::chg_costs,
  //                              MS_dbl_rngd::args() );
  //
  // register_method< OTBlock >( "OTBlock::chg_costs", &OTBlock::chg_costs,
  //                              MS_dbl_sbst::args() );
  //
  // register_method< OTBlock >( "OTBlock::chg_ucaps", &OTBlock::chg_ucaps,
  //                              MS_dbl_rngd::args() );
  //
  // register_method< OTBlock >( "OTBlock::chg_ucaps", &OTBlock::chg_ucaps,
  //                              MS_dbl_sbst::args() );
  //
  // register_method< OTBlock >( "OTBlock::chg_dfcts", &OTBlock::chg_dfcts,
  //                              MS_dbl_rngd::args() );
  //
  // register_method< OTBlock >( "OTBlock::chg_dfcts", &OTBlock::chg_dfcts,
  //                              MS_dbl_sbst::args() );
  //
  // register_method< OTBlock >( "OTBlock::close_arcs",
  //                              &OTBlock::close_arcs,
  //                              MS_rngd::args() );
  //
  // register_method< OTBlock >( "OTBlock::close_arcs",
  //                              &OTBlock::close_arcs,
  //                              MS_sbst::args() );
  //
  // register_method< OTBlock >( "OTBlock::open_arcs", &OTBlock::open_arcs,
  //                              MS_rngd::args() );
  //
  // register_method< OTBlock >( "OTBlock::open_arcs", &OTBlock::open_arcs,
  //                              MS_sbst::args() );


  register_method< OTBlock , MF_dbl_it , Range >( "OTBlock::chg_costs" ,
						   & OTBlock::chg_costs );

  register_method< OTBlock , MF_dbl_it , Subset && , bool >(
   "OTBlock::chg_costs" , & OTBlock::chg_costs );

  register_method< OTBlock , Range >( "OTBlock::close_arcs" ,
				       & OTBlock::close_arcs );

  register_method< OTBlock , Subset && , bool >( "OTBlock::close_arcs" ,
						  & OTBlock::close_arcs );

  register_method< OTBlock , Range >( "OTBlock::open_arcs" ,
				       & OTBlock::open_arcs );

  register_method< OTBlock , Subset && , bool >( "OTBlock::open_arcs" ,
						  & OTBlock::open_arcs );

  }  // end( static_initialization )

/*--------------------------------------------------------------------------*/

 int p2i_x_s( const Variable * var ) const {
  return( std::distance( x.data() ,
			 static_cast< const ColVariable * >( var ) ) );
  }

 int p2i_ub_s( const Constraint * cns ) const {
  return( std::distance( UB.data() ,
			 static_cast< const LB0Constraint * >( cns ) ) );
  }

 int p2i_e_s( const Constraint * cns ) const {
  return( std::distance( E.data() ,
			 static_cast< const FRowConstraint * >( cns ) ) );
  }

 LinearFunction * get_lfo( void ) {
  #ifdef NDEBUG
   return( static_cast< LinearFunction * >( c.get_function() ) );
  #else
   auto lfo = dynamic_cast< LinearFunction * >( c.get_function() );
   assert( lfo );
   return( lfo );
  #endif
  }

 LinearFunction * get_lfc( FRowConstraint * cnsti )
 {
  #ifdef NDEBUG
   return( static_cast< LinearFunction * >( cnsti->get_function() ) );
  #else
   auto lfc = dynamic_cast< LinearFunction * >( cnsti->get_function() );
   assert( lfc );
   return( lfc );
  #endif
  }

 void guts_of_destructor( void );

 void guts_of_add_Modification( p_Mod mod , ChnlName chnl );

 void compute_conditional_bounds( void );

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG

 void CheckAbsVSPhys( void );
 
#endif

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;  // insert OTBlock in the Block factory

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class( OTBlock ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS OTBlockMod -----------------------------*/
/*--------------------------------------------------------------------------*/
/// derived class from Modification for modifications to a OTBlock
/** Derived class from Modification to describe modifications to a OTBlock.
 *  This is actually "sort of abstract", since it does not say exactly what
 *  is changed, this being demanded to derived classes (which do this in
 *  different ways). Note that it is derived from Modification rather than,
 *  say, BlockMod (which has the same structure) because this is a class of
 *  "physical Modification". This means that a OTBlockMod refers to changes
 *  in the "physical representation" of the OTBlock; the corresponding
 *  changes in the "abstract representation" of the OTBlock are dealt with
 *  by means of "abstract Modification", i.e., derived classes from
 *  AModification (as is BlockMod, which is why OTBlockMod is not derived
 *  from BlockMod). */

class OTBlockMod : public Modification
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*---------------------------- PUBLIC TYPES --------------------------------*/
 /// public enum for the types of OTBlockMod
 
 enum OTB_mod_type {
  eChgCost = 0 ,   ///< change the arc costs
  eOpenArc     ,   ///< open arcs
  eCloseArc    ,   ///< close arcs
  eAddArc      ,   ///< add arcs
  eRmvArc          ///< remove arcs
  };

/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 /// constructor: takes the OTBlock and the type

 OTBlockMod( OTBlock * fblock , int type )
  : f_Block( fblock ) , f_type( type ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~OTBlockMod() = default;   ///< destructor, does nothing

/*-------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

 /// returns the [OT]Block to which the OTBlockMod refers

 Block * get_Block( void ) const override  { return( f_Block ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
 /// accessor to the type of modification

 int type( void ) const { return( f_type ); }

/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the OTBlockMod

 void print( std::ostream &output ) const override {
  output << "OTBlockMod[" << this << "]: ";
  switch( f_type ) {
   case( eChgCost ):  output << "change costs "; break;
   case( eOpenArc ):  output << "open arcs "; break;
   case( eCloseArc ): output << "close arcs "; break;
   case( eAddArc ):   output << "add arcs "; break;
   default:           output << "remove arcs ";
   }
  }

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 OTBlock *f_Block;
               ///< pointer to the OTBlock to which the OTBlockMod refers

 int f_type;   ///< type of Modification

/*--------------------------------------------------------------------------*/

 };  // end( class( OTBlockMod ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS OTBlockRngdMod ---------------------------*/
/*--------------------------------------------------------------------------*/
/// derived from OTBlockMod for "ranged" modifications
/** Derived class from OTBlockMod to describe "ranged"
 * modifications to a OTBlock, i.e., modifications that apply to an interval
 * of either arcs or nodes. */

class OTBlockRngdMod : public OTBlockMod
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 /// constructor: takes the OTBlock, the type, and the range

 OTBlockRngdMod( OTBlock * fblock , int type , Block::Range rng )
  : OTBlockMod( fblock , type ) , f_rng( rng ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~OTBlockRngdMod() = default;   ///< destructor, does nothing

/*-------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

 /// accessor to the range

 Block::c_Range & rng( void ) const { return( f_rng ); }
 
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the OTBlockRngdMod

 void print( std::ostream &output ) const override {
  OTBlockMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
  }

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 Block::Range f_rng;     ///< the range

/*--------------------------------------------------------------------------*/

 };  // end( class( OTBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS OTBlockSbstMod ---------------------------*/
/*--------------------------------------------------------------------------*/
/// derived from OTBlockMod for "subset" modifications
/** Derived class from Modification to describe "subset" modifications to a
 *  OTBlock, i.e., modifications that apply to an arbitrary subset of either
 * the arcs or the nodes. */

class OTBlockSbstMod : public OTBlockMod
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*---------------------- CONSTRUCTOR & DESTRUCTOR --------------------------*/

 ///< constructor: takes the OTBlock, the type, and the subset
 /**< Constructor: takes the OTBlock, the type, and the subset. As the the
  * && tells, nms is "consumed" by the constructor and its resources become
  * property of the OTBlockSbstMod object.
  *
  *   NOTE THAT nms IS REQUIRED TO BE ORDERED IN INCREASING SENSE
  *
  * although this is not checked by the class. */

 OTBlockSbstMod( OTBlock * fblock , int type , Block::Subset && nms )
  : OTBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

 virtual ~OTBlockSbstMod() = default;  ///< destructor, does nothing

/*-------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

 /// accessor to the subset

 Block::c_Subset & nms( void ) const { return( f_nms ); }

/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/
 /// print the OTBlockSbstMod

 void print( std::ostream &output ) const override {
  OTBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
  }

/*--------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 Block::Subset f_nms;   ///< the subset

/*--------------------------------------------------------------------------*/

 };  // end( class( OTBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS OTSolution -----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a solution of a OTBlock
/** The OTSolution class, derived from Solution, represents a solution of a
 * OTBlock, i.e.:
 *
 * - an n*m vector of TNumber for the arc transport plan values;
 * 
 * where n and m are the number of atoms of respectively P and Q. */

class OTSolution : public Solution {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*------------------------------- FRIENDS ----------------------------------*/

 friend OTBlock;  ///< make OTBlock friend

/*---------------- CONSTRUCTING AND DESTRUCTING OTSolution ----------------*/

 explicit OTSolution( void ) { }  /// constructor, it has nothing to do

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void deserialize( const netCDF::NcGroup & group ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~OTSolution() = default;  ///< destructor: it is virtual, and empty

/*------------- METHODS DESCRIBING THE BEHAVIOR OF A OTSolution -----------*/

 void read( const Block * block ) override final;

 void write( Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize a OTSolution into a netCDF::NcGroup
 /** Serialize a OTSolution into a netCDF::NcGroup, with the following
  * format:
  *
  * - The dimension "NumAtomsP" containing the number of atoms of the 
  *   probability distribution P. The dimension
  *   is optional.
  * 
  * - The dimension "NumAtomsQ" containing the number of atoms of the 
  *   probability distribution Q. The dimension
  *   is optional.
  *
  * - The dimension "NumTransPlan" containing the number of arcs. The dimension
  *   is optional, if it is not specified then the corresponding variable
  *   "Potentials" is not read (the OTSolution object does not contain any
  *   flow solution).
  *
  * - The variable "FlowSolution", of type double and indexed over the
  *   dimension NumTransPlan. The variable is optional, if it is not specified
  *   then the OTSolution object does not contain any flow solution.
  *
  * - The variable "Potentials", of type double and indexed over the
  *   dimension NumNodes. The variable is optional, if it is not specified
  *   then the OTSolution object does not contain any node potentials. */
 
 void serialize( netCDF::NcGroup & group ) const override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 OTSolution * scale( double factor ) const override final;

 void sum( const Solution * solution , double multiplier ) override final;

 OTSolution * clone( bool empty = false ) const override final;

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream &output ) const override final {
  output << "OTSolution [" << this << "]: " << v_x.size() << " transport plan " 
    << std::endl;
  }

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 OTBlock::Vec_TNumber v_x;   ///< the transport plan

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( OTSolution ) )

/** @} end( group( OTBlock_CLASSES ) ) ------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* OTBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File OTBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/

/*
TODO list

 - Initial admissible transport plan is \pi defined by \pi_{ij} = p_iq_j. Should
    define a helper somewhere?

*/
