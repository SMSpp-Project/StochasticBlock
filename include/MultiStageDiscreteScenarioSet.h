/*--------------------------------------------------------------------------*/
/*------------------ File MultiStageDiscreteScenarioSet.h ------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *concrete* class MultiStageDiscreteScenarioSet, an
 * implementation of MultiStageScenarioGenerator for the case where the
 * multistage scenarios are a *general* (history-dependent) discrete
 * stochastic process, given explicitly as a scenario *tree*: each node at
 * stage t carries a realization x_t of the random variable X_t | H_t and a
 * conditional probability P( x_t | H_t ), and its children are the possible
 * realizations of X_{t+1} | H_{t+1}. This is the discrete counterpart of the
 * "general" multistage case (as opposed to the stage-independent one handled
 * by IndependentMultiStageScenarioGenerator): the realizations of a stage
 * genuinely depend on which branch of the tree one is in, which is exactly
 * what makes the scenario tree branch (rather than recombine into the
 * "fishbone" of a stage-independent process).
 *
 * Besides the legacy single-cursor interface inherited from
 * MultiStageScenarioGenerator (a MultiStageDiscreteScenarioSet *is* a valid
 * MultiStageScenarioGenerator, walkable stage-by-stage by a single consumer),
 * this class also exposes a *view*-based access path: make_view( node )
 * returns a lightweight, read-only ScenarioGenerator pinned to one tree node
 * that iterates the children of that node as an ordinary single-stage pool.
 * Because the tree data is immutable once built, any number of views can be
 * alive at the same time, each with its own cursor, so the scenario tree can
 * be traversed (and the corresponding Block tree built) concurrently, on
 * multiple threads, without contention. A consumer that just wants "the
 * scenarios of this node" (e.g. a TwoStageStochasticBlock) consumes a view
 * as a plain ScenarioGenerator, without being aware that it sits inside a
 * larger tree.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __MultiStageDiscreteScenarioSet
 #define __MultiStageDiscreteScenarioSet
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BlockSolverConfig.h"
#include "ScenarioGenerator.h"

#include <algorithm>
#include <memory>

#include <vector>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*------------- CLASS MultiStageDiscreteScenarioSet ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// concrete MultiStageScenarioGenerator for a general discrete scenario tree
/** MultiStageDiscreteScenarioSet implements the general (history-dependent)
 * discrete multistage case as an explicit scenario *tree*. The process
 * ( X_0 , ... , X_{T-1} ) is described by a rooted tree where:
 *
 *  - the root is the (unique, deterministic) stage-0 node, with empty
 *    history H_0;
 *
 *  - a node at stage t holds a realization x_t of X_t | H_t (a vector of
 *    length d_t) and a conditional probability P( x_t | H_t ); its children
 *    are the possible realizations of X_{t+1} | ( H_t , x_t );
 *
 *  - a root-to-leaf path is a full scenario ( x_0 , ... , x_{T-1} ), whose
 *    probability is the product of the conditional probabilities along it.
 *
 * The conditional probabilities of the children of any node must sum to 1.
 * Stage-independence is the special case where every node at stage t has the
 * same children (same data, same probabilities) regardless of its history;
 * that case is better served by IndependentMultiStageScenarioGenerator, which
 * stores it without the combinatorial blow-up. This class is for when the
 * children genuinely depend on the branch.
 *
 * ### Reading the tree: views
 *
 * The tree is read exclusively through MultiStageScenarioGenerator::View,
 * here implemented by TreeView: a small read-only ScenarioGenerator pinned
 * to one node, which iterates *the children of that node* as an ordinary
 * single-stage pool (with the conditional probabilities as the pool
 * weights), and which moves through the tree with descend() / climb().
 * root_view() returns the one pinned at the root; the generator itself,
 * seen as a plain ScenarioGenerator, behaves exactly like it.
 *
 * The tree data is immutable after deserialize() (and after any global
 * reduction, see below), and each view carries its own cursor, so an
 * arbitrary number of views can be alive and be used *at the same time*.
 * This is what allows building the corresponding Block tree concurrently,
 * one view per subtree, on as many threads as one likes, and what lets a
 * node-local consumer (a TwoStageStochasticBlock) treat its slice of the
 * tree as a plain ScenarioGenerator without knowing about the rest.
 *
 * ### Reducing a node is not selecting scenarios
 *
 * @warning What init_representative_pool( K ) does here is a *local
 *          heuristic*, and a placeholder: it selects K of the children of
 *          one node, judging them by a measure of goodness applied to them
 *          alone. That is not scenario selection. In a multi-stage problem
 *          a scenario is a whole root-to-leaf path, not a node, so asking
 *          for the K best scenarios is a question about paths, whose answer
 *          may well be a tree of a different shape, rather than the same
 *          tree with some branches pruned level by level. Note in
 *          particular that the K asked for here is a number of children of
 *          one node, and not a number of scenarios: the two coincide only
 *          in the last stage of a tree with a single node per earlier one.
 *
 * The two are told apart by where the reduction is asked for.
 * init_representative_pool( K ) called on a View acts on the pool of the
 * node it is pinned at, and is the local heuristic above; called on the
 * generator it acts on the whole object, which is where a global reduction
 * rebuilding the tree belongs, that being the only level holding every
 * path. Today the latter delegates to the View pinned at the root, so the
 * two coincide, but the entry point is there for the real algorithm to fill
 * in, and the generation counters already do the right thing: a reduction
 * at the root invalidates every View below it, which is exactly what has to
 * happen when the shape of the tree changes.
 *
 * The selection itself is not implemented here: the children of the node
 * are handed to a DiscreteScenarioSet, which reduces them through a
 * ScenarioReductionBlock and whatever Solver the BlockSolverConfig passed
 * to set_solver_config() attaches to it. The reduction algorithms thus stay
 * where they belong, in a Solver, and the single-stage and the multi-stage
 * case share them; with no Solver configured the fall-back is the same
 * baseline selection the single-stage case uses.
 *
 * Which is also why what those Solver do is a heuristic here and not the
 * answer: they reduce a flat set of scenario vectors, judging them by a
 * matrix of pairwise distances between the vectors themselves. On the
 * children of one node that is exactly the right question. On a tree it is
 * not, since two paths sharing a prefix cannot be told apart before the
 * prefix ends, so the distance that a multi-stage reduction needs is one
 * that accounts for that, and the shape of the resulting tree is part of
 * its answer rather than an input to it. Feeding the flattened paths to a
 * single-stage algorithm does return something, but it answers a different
 * question.
 *
 * Reducing a node does not throw its children away: they remain its
 * universe, only the pool that the views iterate is restricted, and the
 * conditional probability of each discarded child is accumulated onto the
 * representative it has been assigned to. The reduction is therefore
 * reversible, calling it again with INFScenario restores the whole set of
 * children, and it is local: reducing a node says nothing about the pools
 * of its own children, which the caller reduces by descending into them.
 *
 * ### netCDF format
 *
 * A group containing a MultiStageDiscreteScenarioSet has a mandatory string
 * attribute "type" = "MultiStageDiscreteScenarioSet" plus:
 *
 *  - dimension "NumberStages" (size T > 0);
 *
 *  - dimension "NumberNodes" (size N >= 1, the total number of tree nodes,
 *    root included);
 *
 *  - variable "NodeStage" [NumberNodes]: the stage t of each node;
 *
 *  - variable "NodeParent" [NumberNodes]: the parent node index of each node
 *    (the root, conventionally node 0, has parent = NumberNodes as a
 *    "no parent" marker);
 *
 *  - variable "NodeProbability" [NumberNodes]: P( node | parent ) for each
 *    node (the root has probability 1);
 *
 *  - dimension "ScenarioDataSize" (size D) and variable "NodeData"
 *    [NumberNodes][ScenarioDataSize]: row n holds the realization x_t of
 *    node n, padded to D = max_t d_t (only the first d_{stage(n)} entries
 *    are meaningful). The per-stage sizes d_t are given by variable
 *    "StageScenarioSize" [NumberStages].
 *
 * Nodes are required to be listed parent-before-child (a topological order),
 * so that the tree can be reconstructed in a single pass. */

class MultiStageDiscreteScenarioSet : public MultiStageScenarioGenerator
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

 /// type for indexing a node of the scenario tree
 using NodeIndex = unsigned int;

 /// constexpr marking "no node" (e.g. the parent of the root)
 static constexpr NodeIndex InvalidNode = Inf< NodeIndex >();

/** @} ---------------------------------------------------------------------*/
/*--- CONSTRUCTING AND DESTRUCTING MultiStageDiscreteScenarioSet -----------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing MultiStageDiscreteScenarioSet
 *  @{ */

 MultiStageDiscreteScenarioSet( void )
  : MultiStageScenarioGenerator() { }

/*--------------------------------------------------------------------------*/

 virtual ~MultiStageDiscreteScenarioSet() = default;

/*--------------------------------------------------------------------------*/

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*------------------------- VIEW-BASED ACCESS ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name View-based access to the tree
 *  @{ */

 /// create a View pinned at the root (the MSSG view entry point)
 /** Implements MultiStageScenarioGenerator::root_view(): returns a View (a
  * TreeView) pinned at the root, whose pool are the first-stage
  * realizations. From it the whole tree is walked via View::descend() and
  * View::climb(). This is the only way to read the tree, and it is
  * generator-agnostic: consumers (e.g. MSSB) hold a general
  * MultiStageScenarioGenerator rather than this concrete type. */

 [[nodiscard]] std::unique_ptr< View > root_view( void ) const override;

/** @} ---------------------------------------------------------------------*/
/*------ METHODS INHERITED FROM (Multi-Stage)ScenarioGenerator -------------*/
/*--------------------------------------------------------------------------*/
/** @name Functions inherited from (Multi-Stage)ScenarioGenerator
 *
 * Seen as a plain ScenarioGenerator, a MultiStageDiscreteScenarioSet is its
 * own root view: these methods all read the realizations of the first
 * random variable, i.e. the children of the root. Reaching any other
 * position in the tree is done with root_view() and View::descend().
 *  @{ */

 void set_seed( unsigned long seed = 0 ) override { }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] ScenarioIndex get_support_size( void ) override;

/*--------------------------------------------------------------------------*/

 [[nodiscard]] ScenarioSize get_scenario_size( void ) const override;

/*--------------------------------------------------------------------------*/

 [[nodiscard]] Scenario get_current_scenario( void ) const override;

/*--------------------------------------------------------------------------*/

 [[nodiscard]] double get_current_scenario_probability( void ) const override;

/*--------------------------------------------------------------------------*/

 bool next_scenario( void ) override;

/*--------------------------------------------------------------------------*/

 void reset_pool( void ) override;

/*--------------------------------------------------------------------------*/

 [[nodiscard]] bool is_pool_initialized( void ) const override {
  return( ! f_nodes.empty() );
  }

/*--------------------------------------------------------------------------*/

 [[nodiscard]] StageIndex get_stage_number( void ) override {
  return( f_number_stages );
  }

/*--------------------------------------------------------------------------*/

 // init_random_pool() acts on the pool of the root view;
 // init_representative_pool() reduces the WHOLE tree (see class comments).

 void init_random_pool( ScenarioIndex size = INFScenario ) override;

 void init_representative_pool( ScenarioIndex size = INFScenario ) override;

/*--------------------------------------------------------------------------*/
 /// set the BlockSolverConfig that reduces the pool of a node
 /** Sets the BlockSolverConfig that init_representative_pool() applies to
  * the ScenarioReductionBlock it reduces a node with; ownership of \p
  * solver_config is transferred. With none set, the selection falls back to
  * the baseline one, exactly as in the single-stage case. */

 void set_solver_config( BlockSolverConfig * solver_config ) {
  f_solver_config.reset( solver_config );
  }

/*--------------------------------------------------------------------------*/
 /// set the Block the scenarios of this tree are those of
 /** Stores the Block, which init_representative_pool() passes on to the
  * ScenarioReductionBlock together with the scenarios of the node being
  * reduced, since the cost-aware reduction algorithms need it to evaluate
  * them; ownership is *not* transferred. */

 void set_Block( Block * block ) override {
  f_stochastic_block = block;
  }

/** @} ---------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE TYPES ---------------------------------*/
/*--------------------------------------------------------------------------*/

 /// one node of the scenario tree
 struct Node {
  StageIndex stage;                   ///< stage t of this node
  NodeIndex parent;                   ///< parent node (InvalidNode for root)
  double probability;                 ///< P( this | parent )
  std::vector< double > data;         ///< realization x_t (length d_t)
  std::vector< NodeIndex > children;  ///< child node indices
  };

/*--------------------------------------------------------------------------*/
 /// the pool of a node, i.e. which of its children the views iterate
 /** The tree read from the netCDF file is immutable, while the pool of a
  * node is not: init_representative_pool() restricts it to K of the
  * children, with their probabilities re-distributed, and restores it on
  * request. Keeping it apart from the Node is what makes it possible to
  * reduce a node through a View, which only holds a const pointer to the
  * generator, without pretending that the tree itself changes. */

 struct Pool {
  std::vector< NodeIndex > selected;  ///< the children currently in the pool
  std::vector< double > weights;      ///< their conditional probabilities
  unsigned long generation = 0;       ///< bumped by every (re-)definition
  };

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS TreeView -----------------------------------*/
/*--------------------------------------------------------------------------*/
 /// read-only single-stage view of the children of a tree node
 /** A TreeView is the object returned by
  * MultiStageDiscreteScenarioSet::root_view(). It is a full-fledged
  * (single-stage) ScenarioGenerator whose "pool" is the children of one
  * node, iterated in their natural order with the conditional
  * probabilities as weights, and it moves to another node of the tree with
  * descend() / climb(). It holds a non-owning pointer to its parent
  * MultiStageDiscreteScenarioSet and a private cursor, so distinct
  * TreeView-s never interfere: this is what makes concurrent traversal
  * safe. The only thing a TreeView writes is the pool of the node it is
  * pinned at, through init_representative_pool() (see the class comments of
  * MultiStageDiscreteScenarioSet).
  *
  * Since the tree is materialised upfront, the pool of a TreeView is always
  * initialised and a TreeView is valid for its whole life, save for a global
  * reduction actually rebuilding the tree: a TreeView records the generation
  * of each node of the pinned history, so that it detects a node of it
  * having been regenerated and refuses to read stale data [see
  * MultiStageScenarioGenerator::View::is_valid()]. */

 class TreeView : public MultiStageScenarioGenerator::View
 {
  public:

   TreeView( const MultiStageDiscreteScenarioSet * parent , NodeIndex node )
    : f_parent( parent ) , f_node( node ) , f_pos( 0 ) {
    // record the generation of each node of the pinned history, root first
    for( auto n = node ; n != InvalidNode ; n = parent->get_parent( n ) )
     v_stamps.push_back( parent->get_node_generation( n ) );
    std::reverse( v_stamps.begin() , v_stamps.end() );
    }

   void deserialize( const netCDF::NcGroup & ) override {
    throw( std::logic_error(
     "MultiStageDiscreteScenarioSet::TreeView::deserialize: a view is not "
     "independently (de)serializable" ) );
    }

   void set_seed( unsigned long seed = 0 ) override { }

   [[nodiscard]] ScenarioSize get_scenario_size( void ) const override;

   // number of scenarios in this view's pool = children of the node
   [[nodiscard]] ScenarioIndex get_support_size( void ) override;

   [[nodiscard]] Scenario get_current_scenario( void ) const override;

   [[nodiscard]] double get_current_scenario_probability( void )
    const override;

   bool next_scenario( void ) override;

   void reset_pool( void ) override { f_pos = 0; }

   [[nodiscard]] bool is_pool_initialized( void ) const override {
    return( true );
    }

   // the pool of the node this view is pinned at is the one thing a view
   // writes; random sampling of the children is not supported (yet)
   void init_random_pool( ScenarioIndex = INFScenario ) override { }
   void init_representative_pool( ScenarioIndex size = INFScenario )
    override;

   // history-pinned View interface: move through the tree
   bool descend( void ) override;
   bool climb( void ) override;
   [[nodiscard]] std::unique_ptr< View > clone( void ) const override;
   [[nodiscard]] StageIndex stage( void ) const override;
   [[nodiscard]] bool is_valid( void ) const override;

  private:

   const MultiStageDiscreteScenarioSet * f_parent;  ///< owning tree
   NodeIndex f_node;          ///< node whose children this view iterates
   std::size_t f_pos;         ///< cursor into that node's children

   /// generation of each node of the pinned history, root first
   std::vector< unsigned long > v_stamps;

   /// throw if some node of the pinned history has been regenerated
   void check_valid( const char * method ) const;

   [[nodiscard]] const std::string & private_name( void ) const override;

 };   // end( class TreeView )

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
 // direct access to the node store: the tree is read from the outside only
 // through a View, so these are for the TreeView (and this class) alone.

 /// index of the (unique) root node, i.e. the stage-0 node
 [[nodiscard]] NodeIndex get_root( void ) const { return( 0 ); }

/*--------------------------------------------------------------------------*/
 /// the children of node \p n (the realizations of the next stage there)
 [[nodiscard]] const std::vector< NodeIndex > & get_children(
						   NodeIndex n ) const {
  return( f_nodes.at( n ).children );
  }

/*--------------------------------------------------------------------------*/
 /// the parent of node \p n, InvalidNode if \p n is the root
 [[nodiscard]] NodeIndex get_parent( NodeIndex n ) const {
  return( f_nodes.at( n ).parent );
  }

/*--------------------------------------------------------------------------*/
 /// the stage t in 0, ..., T-1 of node \p n
 [[nodiscard]] StageIndex get_node_stage( NodeIndex n ) const {
  return( f_nodes.at( n ).stage );
  }

/*--------------------------------------------------------------------------*/
 /// the conditional probability P( node \p n | its parent )
 [[nodiscard]] double get_node_probability( NodeIndex n ) const {
  return( f_nodes.at( n ).probability );
  }

/*--------------------------------------------------------------------------*/
 /// how many times the children of node \p n have been redefined
 /** The generation of a node is bumped every time its children are (re-)
  * defined, which is how a TreeView tells that the sub-tree it is pinned
  * into is no longer the one it was created on. In this class the tree is
  * read from a netCDF file and never changes afterwards, so generations
  * only move if a global reduction rebuilds the tree. */

 [[nodiscard]] unsigned long get_node_generation( NodeIndex n ) const {
  return( v_pool.at( n ).generation );
  }

/*--------------------------------------------------------------------------*/
 /// the children of node \p n that are currently in its pool
 [[nodiscard]] const std::vector< NodeIndex > & get_pool(
						   NodeIndex n ) const {
  return( v_pool.at( n ).selected );
  }

/*--------------------------------------------------------------------------*/
 /// the probability of the \p i-th child in the pool of node \p n
 [[nodiscard]] double get_pool_weight( NodeIndex n , std::size_t i ) const {
  return( v_pool.at( n ).weights.at( i ) );
  }

/*--------------------------------------------------------------------------*/
 /// reset the pool of node \p n to all its children, with no reduction
 void reset_node_pool( NodeIndex n ) const;

/*--------------------------------------------------------------------------*/
 /// restrict the pool of node \p n to \p size representative children
 /** Selects \p size representatives among the children of node \p n by
  * handing them to a DiscreteScenarioSet and reducing that through a
  * ScenarioReductionBlock with the configured Solver [see
  * set_solver_config()], then bumps the generation of the node, which
  * invalidates every View pinned in the sub-tree below it. With \p size ==
  * INFScenario the whole set of children is restored. */

 void reduce_node( NodeIndex n , ScenarioIndex size ) const;

/*--------------------------------------------------------------------------*/
 /// the realization x_t stored at node \p n
 [[nodiscard]] Scenario get_node_data( NodeIndex n ) const {
  const auto & nd = f_nodes.at( n );
  return( Scenario( nd.data.data() , nd.data.size() ) );
  }

/*--------------------------------------------------------------------------*/
 /// the root view backing the ScenarioGenerator face of this object
 /** Returns the TreeView pinned at the root that implements the inherited
  * single-stage ScenarioGenerator methods; it is created on first use and
  * reset by deserialize(). Throws if the tree is empty. */

 [[nodiscard]] TreeView & self_view( void ) const;

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE FIELDS --------------------------------*/
/*--------------------------------------------------------------------------*/

 /// the tree nodes, node 0 is the root, listed parent-before-child
 std::vector< Node > f_nodes;

 /// number of stages T
 StageIndex f_number_stages = 0;

 /// per-stage realization size d_t
 std::vector< ScenarioSize > f_stage_size;

 /// the root view implementing the ScenarioGenerator face of this object
 mutable std::unique_ptr< TreeView > f_self_view;

 /// per-node pool state; mutable, unlike the tree the nodes describe
 mutable std::vector< Pool > v_pool;

 /// the BlockSolverConfig that attaches the Solver doing the reduction
 std::unique_ptr< BlockSolverConfig > f_solver_config;

 /// the Block the scenarios refer to, not owned
 Block * f_stochastic_block = nullptr;

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

};   // end( class MultiStageDiscreteScenarioSet )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

#endif  /* MultiStageDiscreteScenarioSet.h included */

/*--------------------------------------------------------------------------*/
/*---------------- End File MultiStageDiscreteScenarioSet.h ----------------*/
/*--------------------------------------------------------------------------*/
