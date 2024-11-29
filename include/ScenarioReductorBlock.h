/*--------------------------------------------------------------------------*/
/*-------------------- File ScenarioReductorBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class ScenarioReductorBlock, which implements
 * the Block concept [see Block.h] for the continuous scenario reduction
 * optimization problem.
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

#ifndef __ScenarioReductor
 #define __ScenarioReductor  /* self-identification: #endif at the file end */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"

#include "QuadFunction.h"

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "Solution.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
 class ScenarioReductorBlock;     // forward declaration ScenarioReductorBlock

 class ScenarioReductorSolution;  // fwd declaration ScenarioReductorSolution


/*--------------------------------------------------------------------------*/
/*---------------- ScenarioReductorBlock-RELATED TYPES ---------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup ScenarioReductorBlock_TYPES ScenarioReductorBlock-related types
 *  @{ */

/** @}  end( group( ScenarioReductorBlock_TYPES ) ) */ 

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup ScenarioReductorBlock_CLASSES Classes in ScenarioReductorBlock.h
 *  @{ */


/** @} end( group( ScenarioReductorBlock_CLASSES ) ) */
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS ScenarioReductorBlock -----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for continuous scenario reduction
/** The ScenarioReductorBlock class implements the Block concept [see Block.h] 
 * for the continuous scenario reduction problem.
 * 
 * In a nutshell, the goal is to approximate an input discrete probability 
 * distribution P with n atoms by a discrete probability distribution Q with a 
 * given support size m << n.
 * 
 * We denote by (p_i)_{i \in I} (resp. (q_j)_{j \in J}) the probability weights
 * of the atoms (x_i)_{i \in I} of P (resp. (y_j)_{j \in J} of Q). The atoms x_i
 * and y_j live in the same normed space, we write ||.|| its norm. 
 * 
 * Informally, given a budget m > 0, we want to find the best probability 
 * distribution Q with m atoms that minimizes the k-Wasserstein distance between
 * P and Q.
 * 
 * Formally, given k >= 1, writing c[ i, j ] = || x_i - y_j ||^k we aim to 
 * solve the so-called Continuous Scenario Reduction (CSR) optimization problem
 * \f[
 *  \min_{\pi, (y_j), (q_j)} \sum_{ (i, j) \in I \times J } c[ i, j ] \pi[ i, j ]
 * \f]
 * \f[
 *  \sum_{ j \in J } \pi[ i, j ] = p_i \quad i \in I              (1)
 *  \sum_{ i \in I } \pi[ i, j ] = q_j \quad j \in J              (2)
 *  \pi[ i , j ] \geq 0,
 * \f]
 * where \pi \in \mathbb{R}^{I \times J} is called the transport plan. The 
 * coefficient (i, j) of the transport plan encodes the proportion of the i-th
 * atom of I that is sent to the j-th atom of J and vice versa. Thus, the 
 * transport plan \pi has to satisfy the mass conservation constraints (1) and 
 * (2). 
 * 
 * We reformulate the above optimization problem as a non-convex QP,
 * \f[
 *  \min_{\pi, y_j, q_j, c_{ij}} \sum_{ (i, j) \in I \times J } c_{ij} \pi[i,j]
 * \f]
 * \f[
 *  \sum_{ j \in J } \pi[ i, j ] = p_i \quad i \in I              
 *  \sum_{ i \in I } \pi[ i, j ] = q_j \quad j \in J              
 *  \pi[ i , j ] \geq 0
 *  c_{ij} \geq || x_i - y_j ||^k \quad (i, j) \in I \times J.    (3)
 * \f]
 * The above reformulation with the additional epigraphical variables c_{ij} is 
 * the one encoded by ScenarioReductorBlock as the abstract representation of 
 * ScenarioReductorBlock's optimization problem.
 * 
 * TODO:
 * When the support (y_j) of the probability distribution Q is constrained to
 * be a subset of the support (x_i) of the probability distribution P, then
 * CSR is usually called Discrete Scenario Reduction (DSR). We (=Romain?) 
 * provide an "intermediary" solver implementing Dupacova and al. algorithm 
 * which aims to solve DSR, see [DupacovaDSRSolver.h]. It extracts the physical
 * data out of ScenarioReductorBlock and solves the associated DSR instead of 
 * the full CSR optimization problem.
 * 
 * Note that the DSR solution can be seen an intermediary step to solving the
 * CSR optimization problem.
 * 
 * TODO:
 * For now only the case of the euclidean norm is considered. Should extend
 * ScenarioReductorBLock to arbitrary cost functions. Important use cases for
 * such extension are the Fortet-Mourier metric and the problem-dependant cost
 * of Bertsimas and Mundru. It should be the matter of changing the Function
 * in associated with the FRowConstraint (3).
 * Maybe have one more physical internal variable that can be loaded. And when
 * the FRowConstraint is constructed use that additional variable to change
 * the associated Function? 
 * 
 * Should ask Antonio on how to design this I think.
 * 
 * TODO:
 * generate_abstract_variables()
 * guts_of_destructor()
 * */

class ScenarioReductorBlock : public Block
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:
/** @name Public types 
 @{ */

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

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

/// type to hold the atoms of a probability distribution
typedef double ANumber; 
typedef const ANumber c_ANumber;

typedef boost::multi_array< ANumber, 2 > Mat_ANumber;
typedef const Mat_ANumber c_Mat_ANumber;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/// type to hold the weights of each atom of a probability distribution
typedef double WNumber;
typedef const WNumber c_WNumber;

typedef std::vector< WNumber > Vec_WNumber;
typedef const Vec_WNumber c_Vec_WNumber;

/** @} ---------------------------------------------------------------------*/
/*------------------------------- FRIENDS ----------------------------------*/
/*--------------------------------------------------------------------------*/

// currently no friend (haha), maybe ScenarioReductorSolution if needed?

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of ScenarioReductorBlock, taking a pointer to father Block
 /** Constructor of ScenarioReductorBlock. It accepts a pointer to the father 
  * Block, which can be of any type, defaulting to nullptr so that this can 
  * also be used as the void constructor. */

explicit ScenarioReductorBlock( Block *father = nullptr ) 
  : Block{father}, AR{0}, N{0} , M{0}, k{0} {}

/// destructor of ScenarioReductorBlock: deletes the abstract representation

virtual ~ScenarioReductorBlock() { guts_of_destructor(); }

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * 
 *  *  @{ */

/// loads the CSR instance from memory by copying
/** Loads the CSR instance from memory by copying the data from the input.
* 
*  The parameters are what you expect:
*
* - N   number of atoms in probability distribution P
* 
* - M   number of desired atoms in probability distribution Q
*
* - dim number of dimensions of the euclidean space in which the atoms live
*
* - k   order of the considered k-Wasserstein distance between P and Q
*
* - vP  vector of the atoms of the probability distribution P, size must be n
*
* - wP  vector of weights of the probability distribution P, size must be n */

 void load( Index n,
            Index m,
            Index dim,
            c_Mat_ANumber & atoms,
            c_Vec_WNumber & weights,
            Index k = 2 );

/*--------------------------------------------------------------------------*/
/// loads the CSR instance from memory by moving
/** Loads the CSR instance from memory by moving the data from the input.
* 
*  The parameters are what you expect:
*
* - N   number of atoms in probability distribution P
* 
* - M   number of desired atoms in probability distribution Q
*
* - dim number of dimensions of the euclidean space in which the atoms live
*
* - k   order of the considered k-Wasserstein distance between P and Q
*
* - vP  vector of the atoms of the probability distribution P
*
* - wP  vector of probability weights of the probability distribution P */

 void load( Index n,
            Index m,
            Index dim,
            c_Mat_ANumber && atoms,
            c_Vec_WNumber && weights,
            Index k = 2 );

/*--------------------------------------------------------------------------*/
 /// load instance from .txt file  
 /** Loads a ScenarioReductorBlock out of std::istream. The format is the
  * following, with each element being separated by whitespaces and possibly
  * comments:
  *
  * - N number of atoms in P
  *
  * - M number of desired atoms of Q
  * 
  * - dim number of dimensions of the euclidean space in which the atoms live
  *
  * - for i = 1 to N; for d = 1 to dim; d-th component of the i-th atom of P
  *
  * - for i = 1 to N; weight of atom i of P
  *
  * Since there is only one supported input format, \p frmt is ignored.
  *
  * Like load( memory ), if there is any Solver attached to this 
  * ScenarioReductorBlock then a NBModification (the "nuclear option") is 
  * issued. */

 void load( std::istream & input , char frmt = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// extends Block::deserialize( netCDF::NcGroup )
 /** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
  * a ScenarioReductorBlock. Besides what is managed by serialize() method of 
  * the base Block class, the group should contains the following:
  *
  * - the dimension "N" containing the number of atoms of the probability 
  *   distribution P
  *
  * - the dimension "M" containing the desired number of atoms of the 
  *   probability distribution Q
  * 
  * - the dimension "dim" containing the dimension of the atoms' associated
  *   euclidean space
  *
  * - the variable ""
  *
  * - the variable "Profits", of type double and indexed over the dimension
  *   "NItems"; the i-th entry of the variable is assumed to contain the 
  *   profit of the i-th item;
  * 
  * All dimensions and variables are mandatory. */

 void deserialize( const netCDF::NcGroup & group ) override;

/** @} ---------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE DATA OF ScenarioReductorBlock ---------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of ScenarioReductorBlock
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// get the number of atoms of the probability distribution P

 Index get_N( void ) const { return( N ); }

/*--------------------------------------------------------------------------*/
 /// get the dimension of the atoms' space

 Index get_dim( void ) const { return( dim ); }


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

// "Physical" fields
 Index N;               ///< the number of atoms of P
 Index M;               ///< the maximum number of atoms of Q
 Index k;               ///< the order of the Wasserstein distance
 Index dim;             ///< dimension of the atoms' space
 Mat_ANumber atomsP;      ///< vector holding the atoms of P
 Vec_WNumber weightsP;  ///< vector holding the weights of P

// Abstract representation related fields
 unsigned char AR;               ///< bit-wise coded: what abstract is there

 static constexpr unsigned char HasVar = 1;
 ///< first bit of AR == 1 if the Variable have been constructed
 static constexpr unsigned char HasObj = 2;
 ///< second bit of AR == 1 if the Objective has been constructed
 static constexpr unsigned char HasMas = 4;
 ///< third bit of AR == 1 if the mass conservations have been constructed
 static constexpr unsigned char HasEpi = 8;
 ///< fourth bit of AR == 1 if epigraphical constraints have been constructed

 std::vector< ColVariable > TP;         ///< TP variables
 std::vector< FRowConstraint > massP;   ///< TP mass constraints on P
 std::vector< FRowConstraint > massQ;   ///< TP mass constraints on Q
 
 std::vector< ColVariable > Cost;       ///< Cost variables
 std::vector< FRowConstraint > epiCost; ///< Cost epigraphical constraints  

 FRealObjective objective;               ///< the objective function

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void guts_of_destructor( void );

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;  // insert ScenarioReductorBlock in the factory

/*--------------------------------------------------------------------------*/

};  // end( class( ScenarioReductorBlock ) )

/*--------------------------------------------------------------------------*/

};  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* ScenarioReductor.h included */

/*--------------------------------------------------------------------------*/
/*----------------- End File ScenarioReductorBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
