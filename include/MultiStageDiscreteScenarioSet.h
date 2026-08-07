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

#include "ScenarioGenerator.h"

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
 * ### Scenario selection is global
 *
 * For a scenario *tree*, scenario reduction is inherently a *global*
 * operation: reducing each stage independently would lose the joint
 * structure of the tree. Hence init_representative_pool( K ) on the
 * MultiStageDiscreteScenarioSet itself reduces the *whole tree* at once
 * (rebuilding an equivalent, smaller immutable tree), unlike the
 * per-position contract documented on the scalar forms in
 * ScenarioGenerator.h. Consistently, TreeView::init_*_pool() is a no-op: a
 * view never reduces on its own, it only reads the (already globally
 * reduced) tree. The intended pattern is therefore: build / deserialize the
 * tree, optionally call init_representative_pool( K ) once on the whole
 * generator, then fan out views.
 *
 * @note In this prototype init_representative_pool() performs the *identity*
 *       reduction (the tree is left as is); the global tree-reduction
 *       algorithm is a separate piece of work. The interface is the point
 *       here: where reduction lives (global, on the generator) and where it
 *       does not (per-view).
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
  * safe. A TreeView reads only; init_*_pool() are no-ops because any
  * scenario reduction is performed globally on the parent (see the class
  * comments of MultiStageDiscreteScenarioSet). */

 class TreeView : public MultiStageScenarioGenerator::View
 {
  public:

   TreeView( const MultiStageDiscreteScenarioSet * parent , NodeIndex node )
    : f_parent( parent ) , f_node( node ) , f_pos( 0 ) { }

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

   // reduction is global, done on the parent: a view never reduces
   void init_random_pool( ScenarioIndex = INFScenario ) override { }
   void init_representative_pool( ScenarioIndex = INFScenario ) override { }

   // history-pinned View interface: move through the tree
   bool descend( void ) override;
   bool climb( void ) override;
   [[nodiscard]] std::unique_ptr< View > clone( void ) const override;
   [[nodiscard]] StageIndex stage( void ) const override;

  private:

   const MultiStageDiscreteScenarioSet * f_parent;  ///< owning tree
   NodeIndex f_node;          ///< node whose children this view iterates
   std::size_t f_pos;         ///< cursor into that node's children

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
