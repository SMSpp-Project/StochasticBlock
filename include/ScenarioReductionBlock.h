/*--------------------------------------------------------------------------*/
/*------------------- File ScenarioReductionBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the ScenarioReductionBlock class, which is a generic Block
 * used to perform scenario reduction. It holds a pointer to a
 * ScenarioGenerator and an inner ScenarioReductionBlockSolution. Any
 * specialised Solver can dynamic_cast<> to ScenarioReductionBlock and write
 * its solution into it. DiscreteScenarioSet then reads the solution and
 * acts upon it.
 *
 * \author Minh Duc Pham \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Minh Duc Pham
 */
/*--------------------------------------------------------------------------*/

#ifndef __ScenarioReductionBlock
#define __ScenarioReductionBlock

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"
#include "ScenarioGenerator.h"

#include <vector>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*----------------- CLASS ScenarioReductionBlockSolution -------------------*/
/*--------------------------------------------------------------------------*/
/** Solution produced by a ScenarioReductionSolver and stored inside a
 * ScenarioReductionBlock. It contains the indices of the selected
 * representative scenarios, the assignment of every scenario to its
 * representative, and the aggregated probability weights*/

class ScenarioReductionBlockSolution
{
public:

 using ScenarioIndex = ScenarioGenerator::ScenarioIndex;

 ScenarioReductionBlockSolution() = default;

 /** Indices (in the original set) of the selected representative scenarios*/
 std::vector< ScenarioIndex > selected_indices;

 /** For each scenario i, assignments[i] is the index of the representative
  * it is assigned to. Must have size equal to the total number of scenarios
  * in the ScenarioGenerator*/
 std::vector< ScenarioIndex > assignments;

 /** Aggregated probability weights for each representative scenario.
  * The weights sum to 1.0 and have the same size as selected_indices*/
 std::vector< double > weights;

 /** Returns true if the solution has been written by a Solver*/
 [[nodiscard]] bool is_set() const {
  return( ! selected_indices.empty() );
  }

 /** Clears all solution data*/
 void clear() {
  selected_indices.clear();
  assignments.clear();
  weights.clear();
  }

 };  // end( class ScenarioReductionBlockSolution )

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS ScenarioReductionBlock -----------------------*/
/*--------------------------------------------------------------------------*/
/** A generic Block for scenario reduction. It holds a pointer to a
 * ScenarioGenerator (typically a DiscreteScenarioSet) and an inner
 * ScenarioReductionBlockSolution.
 *
 * Usage pattern:
 *  1. DiscreteScenarioSet creates a ScenarioReductionBlock and passes itself
 *     as the ScenarioGenerator.
 *  2. If a TwoStageStochasticBlock is available (needed by solvers such as
 *     CSSCScenarioReductionSolver), it is set via set_stochastic_block() and
 *     added as the first sub-Block.
 *  3. A specialised Solver is attached and compute() is called.
 *  4. The Solver dynamic_cast<>s to ScenarioReductionBlock, writes its
 *     solution into the ScenarioReductionBlockSolution, and finishes.
 *  5. DiscreteScenarioSet reads the solution and updates its pool. */

class ScenarioReductionBlock : public Block
{

public:

/*--------------------------------------------------------------------------*/
/*------------- CONSTRUCTING AND DESTRUCTING ScenarioReductionBlock --------*/
/*--------------------------------------------------------------------------*/

 /** Constructor
  * @param father  Optional father Block */
 explicit ScenarioReductionBlock( Block * father = nullptr )
  : Block( father ) , f_scenario_generator( nullptr ) {}

 /** Destructor. Does not own the ScenarioGenerator */
 ~ScenarioReductionBlock() override = default;

/*--------------------------------------------------------------------------*/
/*-------------------------- SETTERS ---------------------------------------*/
/*--------------------------------------------------------------------------*/

 /** Set the ScenarioGenerator that provides the scenario data.
  * Ownership is NOT transferred; the caller retains ownership.
  * @param sg  Pointer to the ScenarioGenerator*/
 void set_scenario_generator( ScenarioGenerator * sg ) {
  f_scenario_generator = sg;
  }

 /** Set the TwoStageStochasticBlock used by CSSCScenarioReductionSolver
  * to access the N sub-Blocks and the AbstractPath identifying here-and-now
  * variables. Added as the first sub-Block. Ownership not transferred
  * @param tssb  Pointer to the TwoStageStochasticBlock*/
 void set_stochastic_block( Block * tssb );

 /** Set the StochasticBlock used as a scenario applicator by
  * CSSCScenarioReductionSolver to inject scenario j into sub-Block i
  * during Step 1 of the CSSC algorithm.
  *
  * The applicator wraps the same inner Block as the sub-Blocks of the
  * TwoStageStochasticBlock and carries the DataMappings that know how to
  * map a scenario vector to the Block's parameters (e.g. customer demands
  * for CFL). This avoids any dependency on the specific Block type inside
  * the solver.
  *
  * Ownership is not transferred
  * @param stoch  Pointer to the StochasticBlock applicator*/
 void set_scenario_applicator( Block * stoch ) {
  f_stoch_applicator = stoch;
  }

 /** Set the sub-problem block (e.g. synthetic N×N CFLB) that the solver
  * uses to read nb_atoms, nb_reduced, weights, and distances. */
 void set_sub_problem_block( Block * b ) { f_sub_problem_block = b; }
 [[nodiscard]] Block * get_sub_problem_block() const {
  return( f_sub_problem_block );
  }

 /** Write a solution into the inner ScenarioReductionBlockSolution.
  * Called by a specialised Solver after compute()
  * @param sol  The solution to store. */
 void set_solution( const ScenarioReductionBlockSolution & sol ) {
  f_solution = sol;
  }

/*--------------------------------------------------------------------------*/
/*--------------------------- GETTERS --------------------------------------*/
/*--------------------------------------------------------------------------*/

 /** Returns the ScenarioGenerator. May be nullptr if not set. */
 [[nodiscard]] ScenarioGenerator * get_scenario_generator() const {
  return( f_scenario_generator );
  }

 /** Returns a const reference to the inner solution.
  * Check solution.is_set() before using. */
 [[nodiscard]] const ScenarioReductionBlockSolution & get_solution() const {
  return( f_solution );
  }

 /** Returns a pointer to the TwoStageStochasticBlock sub-Block, or nullptr
  * if it has not been set. */
 [[nodiscard]] Block * get_stochastic_block() const {
  return( v_Block.empty() ? nullptr : v_Block.front() );
  }

 /** Returns a pointer to the StochasticBlock applicator, or nullptr if not
  * set. Used by CSSCScenarioReductionSolver to inject scenario j into
  * sub-Block i during Step 1. */
 [[nodiscard]] Block * get_scenario_applicator() const {
  return( f_stoch_applicator );
  }

/*--------------------------------------------------------------------------*/
/*----------- METHODS FOR PRINTING AND SAVING ------------------------------*/
/*--------------------------------------------------------------------------*/

 void load( std::istream & input , char frmt = 0 ) override {
  throw std::logic_error(
   "ScenarioReductionBlock::load: not implemented" );
  }

 void print( std::ostream & output , char vlvl = 0 ) const override;

 void serialize( netCDF::NcGroup & group ) const override;

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED FIELDS -----------------------------*/
/*--------------------------------------------------------------------------*/

 /** Pointer to the ScenarioGenerator. Not owned by this Block. */
 ScenarioGenerator * f_scenario_generator;

 /** The solution written by a specialised Solver. */
 ScenarioReductionBlockSolution f_solution;

 /** Pointer to the StochasticBlock used as scenario applicator by
  * CSSCScenarioReductionSolver. Not owned. */
 Block * f_stoch_applicator   = nullptr;
 Block * f_sub_problem_block = nullptr;

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

private:

 SMSpp_insert_in_factory_h;

 };  // end( class ScenarioReductionBlock )

}  // end( namespace SMSpp_di_unipi_it )

#endif  /* ScenarioReductionBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------- End File ScenarioReductionBlock.h --------------------*/
/*--------------------------------------------------------------------------*/