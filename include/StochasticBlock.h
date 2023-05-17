/*--------------------------------------------------------------------------*/
/*---------------------- File StochasticBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for the StochasticBlock class, which derives from Block and is
 * meant to turn any Block into its stochastic version.
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Rafael Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __StochasticBlock
#define __StochasticBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"
#include "DataMapping.h"

#include <Eigen/Dense>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup StochasticBlock_CLASSES Classes in StochasticBlock.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS StochasticBlock ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a Block for turning any Block into its stochastic version
/** The StochasticBlock class is meant to represent a Block whose data may be
 * stochastic. The idea is that any Block could have its stochastic version
 * without any changes in its implementation. The StochasticBlock comes to
 * facilitate this process. A StochasticBlock is characterized by the
 * following:
 *
 * 1) It has a pointer to a Block (called here as the inner Block), which is
 *    the Block that is becoming stochastic. This can be any :Block.
 *
 * 2) It is aware that some data of the inner Block is stochastic (and may be
 *    subject to changes) and that the value for this data is represented by a
 *    vector of double. An instance of this vector is called a scenario for
 *    the stochastic data.
 *
 * 3) It has a set of DataMapping. This is used to both identify the
 *    stochastic data in the inner Block and modify the values of this
 *    data. The inner Block may have different types of stochastic data,
 *    located in different sub-Blocks. A DataMapping is meant to represent one
 *    piece of the stochastic data.
 *
 *    To understand how it could be used, we present the following simple
 *    example. Consider a Block B which would be the inner Block of a
 *    StochasticBlock. This Block B may have all sorts of data. We could have
 *    a stochastic version of B by selecting some of its data to become
 *    uncertainties. Let us suppose that B has a vector that represents some
 *    demand over time (from time 0 to T-1). One could establish that the
 *    demand at times 0, 3, and 8 is stochastic. Suppose also that B has a
 *    sub-Block B1 that has a vector representing costs to produce certain
 *    goods. One could also establish that the cost to produce good number 6
 *    is stochastic. In this example, the stochastic data is formed by the
 *    demand of Block B at times 0, 3, and 8, and the cost to produce number 6
 *    of Block B1. DataMapping can then be used to identify each of those
 *    uncertainties.
 *
 *    We could say that a vector of length 4 is a scenario for this stochastic
 *    data such that the first three elements of this vector are values for
 *    the demands and the fourth element is a value for the cost of producing
 *    good number 6. We could use a DataMapping to identify that demand as
 *    stochastic data. To this end, we need two sets. One set identifies what
 *    part of the scenario (indices of that vector) is associated with the
 *    demand and the second set identifies what part of the demand of B is
 *    stochastic. In this example, the first set would be S1 = {0, 1, 2}, which
 *    states that the values for the demand are present at positions 0, 1, and
 *    2 of the scenario vector, and the second set would be S2 = {0, 3, 8}, which
 *    states that the demands that are stochastic are those at positions 0, 3,
 *    and 8. We also need a way of modifying this data in B. For this, it
 *    suffices a method in B that can be used to update its demand (see
 *    SimpleDataMapping), let us say one called set_demand( S2 , vector ).
 *
 *    Another DataMapping could be used to identify the cost of producing good
 *    number 6 as stochastic. This DataMapping would have S1 = {4} as the
 *    first set (stating that the value for this data is at position 4 of the
 *    scenario vector) and S2 = {5} as the second one (indicating which good
 *    cost is stochastic). Again, we need a way to update this cost in B1,
 *    which could be done, for instance, if B1 provides a method like
 *    set_cost( S2 , vector ).
 *
 * 4) It provides a method called set_data() which has a scenario (a vector)
 *    as parameter and that sets the data of the inner Block according to the
 *    set of DataMapping present in the StochasticBlock. In the example above,
 *    when the set_data() method of StochasticBlock is called, the values for
 *    the demands and the good cost are updated with the given values in the
 *    scenario.
 *
 *    In fact, the data of the inner Block does not need to be stochastic in
 *    any sense. What this class provides is a means to set the value of the
 *    data of its inner Block.
 *
 * A StochasticBlock should have a probability distribution (or some kind of
 * partial stochastic process) that describes the uncertainty in it. However,
 * for the moment, it is not supported by this class and this feature will be
 * implemented later on. Typically, an object of this class would be used in
 * conjunction with a scenario generator and the set_data() method of this
 * object would be called to consider a particular scenario.
 */

class StochasticBlock : public Block
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*------------- CONSTRUCTING AND DESTRUCTING StochasticBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing StochasticBlock
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// constructor
 /** Constructs a StochasticBlock having the given \p father Block and \p
  * inner_block inner Block. Both input parameters have a default nullptr
  * value, so that this can be used as the void constructor.
  *
  * @param father The father Block of this StochasticBlock.
  *
  * @param inner_block A pointer to the Block that is becoming stochastic.
  */
 StochasticBlock( Block * father = nullptr ,
                  Block * inner_block = nullptr ) : Block( father ) {
  if( inner_block )
   v_Block.push_back( inner_block );
 }

/*--------------------------------------------------------------------------*/
 /// load the StochasticBlock out of an istream - not implemented yet

 void load( std::istream & input , char frmt = 0 ) override {}

/*--------------------------------------------------------------------------*/
 /// de-serialize a StochasticBlock out of netCDF::NcGroup
 /** The method takes a netCDF::NcGroup supposedly containing all the
  * information required to de-serialize the StochasticBlock, in the format
  * explained in the comments of the serialize() function.
  *
  * @param group a netCDF::NcGroup holding the data in the format described
  *        in the comments to serialize(),
  */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
 /// destructor of StochasticBlock
 /** Destructor of StochasticBlock. It destroys the inner Block (if any),
  * releasing its memory. If the inner Block should not be destroyed then,
  * before this StochasticBlock is destroyed, the pointer to the inner Block
  * must be set to \c nullptr. This can be done by invoking set_inner_block(),
  * passing \c nullptr as a pointer to the new inner Block and \c false to the
  * \c destroy_previous_block parameter. */

 virtual ~StochasticBlock() {
  if( ! v_Block.empty() ) {
   assert( v_Block.size() == 1 );
   delete v_Block.front();
  }
  v_Block.clear();
 }

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// set the (only) sub-Block of this StochasticBlock
 /** This method sets the only sub-Block of this StochasticBlock.
  *
  * @param block the pointer to the Block that is becoming stochastic.
  *
  * @param destroy_previous_block indicates whether the previous inner Block
  *        must be destroyed. The default value of this parameter is \c true,
  *        which means that the previous inner Block (if any) is destroyed and
  *        its allocated memory is released.
  */
 void set_inner_block( Block * block , bool destroy_previous_block = true ) {
  if( ( ! v_Block.empty() ) && ( block == v_Block.front() ) &&
      ( ! destroy_previous_block ) )
   return; // the given Block is already here; silently return

  if( destroy_previous_block && ( ! v_Block.empty() ) ) {
   assert( v_Block.size() == 1 );
   delete v_Block.front();
  }

  v_Block.clear();
  v_Block.push_back( block );

  if( block )
   block->set_f_Block( this );

 if( anyone_there() )
  add_Modification( std::make_shared< NBModification >( this ) );
 }

/*--------------------------------------------------------------------------*/

 /// set the vector of pointers to SimpleDataMappingBase
 /** This method sets the vector of pointers to SimpleDataMappingBase of this
  * StochasticBlock.
  *
  * @param data_mappings the vector of pointers to SimpleDataMappingBase.
  *
  * @param destroy_previous_block indicates whether the previous inner Block
  *        must be destroyed. The default value of this parameter is \c true,
  *        which means that the previous inner Block (if any) is destroyed and
  *        its allocated memory is released.
  */
 void set_data_mappings( std::vector< std::unique_ptr< SimpleDataMappingBase > >
                         && data_mappings ) {
  this->data_mappings = std::move( data_mappings );
 }

/** @} ---------------------------------------------------------------------*/
/*-------------------- Methods for handling Modification -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Modification
 *  @{ */

 void add_Modification( sp_Mod mod , Observer::ChnlName chnl = 0 ) override;

/** @} ---------------------------------------------------------------------*/
/*------------ METHODS FOR Saving THE DATA OF THE StochasticBlock ----------*/
/*--------------------------------------------------------------------------*/
/** @name Saving the data of the StochasticBlock
 *  @{ */

 /// prints the StochasticBlock onto an ostream

 void print( std::ostream & output , char vlvl = 0 ) const override;

/*--------------------------------------------------------------------------*/
 /// serialize a StochasticBlock into a netCDF::NcGroup
 /** Serialize a StochasticBlock into a netCDF::NcGroup, with the following
  * format:
  *
  * - The group "Block", containing the description of the inner Block. This
  *   group is optional. If it is not provided, then the inner Block must be
  *   provided by other means.
  *
  * - The description of a vector of SimpleDataMappingBase, containing the
  *   DataMappings associated with this StochasticBlock. This is optional. For
  *   the time being, each provided DataMapping is expected to be a
  *   SimpleDataMappingBase. Moreover, the inner Block of this StochasticBlock
  *   will serve as the reference Block for both the serialization and
  *   deserialization of each SimpleDataMappingBase.
  */

 virtual void serialize( netCDF::NcGroup & group ) const override;

/**@} ----------------------------------------------------------------------*/
/*-------------- METHODS FOR MODIFYING THE StochasticBlock -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the StochasticBlock
  *  @{ */

/*--------------------------------------------------------------------------*/
 /// sets the (possibly stochastic) data of this StochasticBlock
 /** This function sets the value of the (possibly stochastic) data of this
  * StochasticBlock.
  *
  * @param data The vector containing the values of the data to be set.
  *
  * @param issuePMod Decides if and how a "physical Modification" is issued,
  *        as described in Observer::make_par().
  *
  * @param issueAMod Decides if and how an "abstract Modification" is issued,
  *        as described in Observer::make_par().
  */
 void set_data( const std::vector< double > & data ,
                c_ModParam issuePMod = eNoBlck ,
                c_ModParam issueAMod = eNoBlck ) {
  for( size_t i = 0 ; i < data_mappings.size() ; ++i )
   data_mappings[ i ]->set_data( data.begin() , issuePMod , issueAMod );
 }

/*--------------------------------------------------------------------------*/

 /// sets the (possibly stochastic) data of this StochasticBlock
 /** This function sets the value of the (possibly stochastic) data of this
  * StochasticBlock.
  *
  * @param data An iterator to the first element of the data.
  *
  * @param issuePMod Decides if and how a "physical Modification" is issued,
  *        as described in Observer::make_par().
  *
  * @param issueAMod Decides if and how an "abstract Modification" is issued,
  *        as described in Observer::make_par().
  */
 template< class Iterator >
 void set_data( Iterator data , c_ModParam issuePMod = eNoBlck ,
                c_ModParam issueAMod = eNoBlck ) {
  for( size_t i = 0 ; i < data_mappings.size() ; ++i )
   data_mappings[ i ]->set_data( data , issuePMod , issueAMod );
 }

/*--------------------------------------------------------------------------*/

 /// adds a new SimpleDataMappingBase to this StochasticBlock
 /** This function adds a new SimpleDataMappingBase to the set of
  * SimpleDataMappingBase of this StochasticBlock.
  *
  * @param data_mapping The SimpleDataMappingBase to be added.
  */
 void add_data_mapping( std::unique_ptr< SimpleDataMappingBase > data_mapping ) {
  data_mappings.push_back( std::move( data_mapping ) );
 }

/**@} ----------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE DATA OF THE StochasticBlock -----------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the StochasticBlock
    @{ */

/*--------------------------------------------------------------------------*/

 /// returns the sense of the Objective of this StochasticBlock
 /** This function returns the sense of the Objective of this StochasticBlock,
  * which is defined to be the sense of the Objective of its inner Block. If
  * this StochasticBlock has no inner Block, this function returns
  * Objective::eUndef.
  *
  * @return the sense of the Objective of the inner Block of this
  *         StochasticBlock. */

 int get_objective_sense() const override {
  if( auto inner_block = get_inner_block() )
   return inner_block->get_objective_sense();
  return Objective::eUndef;
 }

/*--------------------------------------------------------------------------*/

 /// returns the vector of pointers to SimpleDataMappingBase
 /** This function returns the vector of pointers to SimpleDataMappingBase
  * that characterizes the data that can be modified through a call to
  * set_data().
  *
  * @return The vector of pointers to SimpleDataMappingBase.
  */

 const std::vector< std::unique_ptr< SimpleDataMappingBase > > &
 get_data_mappings() const {
  return data_mappings;
 }

/*--------------------------------------------------------------------------*/

 /// returns a pointer to the inner Block
 /** This function returns a pointer to the inner Block of this
  * StochasticBlock.
  *
  * @return A pointer to the inner Block of this StochasticBlock.
  */

 Block * get_inner_block() const {
  return v_Block.empty() ? nullptr : v_Block.front();
 }

/**@} ----------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS  ---------------------------*/
/*--------------------------------------------------------------------------*/

 /// The vector of data mappings
 std::vector< std::unique_ptr< SimpleDataMappingBase > > data_mappings;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE FIELDS ------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

};   // end( class StochasticBlock )

/** @} end( group( StochasticBlock_CLASSES ) ) */

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* StochasticBlock.h included */

/*--------------------------------------------------------------------------*/
/*--------------------- End File StochasticBlock.h -------------------------*/
/*--------------------------------------------------------------------------*/
