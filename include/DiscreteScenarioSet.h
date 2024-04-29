/*--------------------------------------------------------------------------*/
/*-------------------- File DiscreteScenarioSet.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file 
 * Header file for the *concrete* class DiscreteScenarioSet that is an 
 * implementation of ScenarioGenerator suited to the case where the input 
 * distribution is contained in a netCDF file as a collection of vectors.
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

#ifndef __DiscreteScenarioSet
    #define __DiscreteScenarioSet
/* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ScenarioGenerator.h"

#include <vector>
#include <random>
#include <span>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*--------------------- CLASS DiscreteScenarioSet --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/

/// DiscreteScenarioSet to sample from a collection of scenarios
/** DiscreteScenarioSet is an implementation of the ScenarioGenerator class.
 * As such, it gives methods to sample from an input distribution and manipulate
 * a scenarioPool. 
 * 
 * In the specific context of DiscreteScenarioSet, the distribution to sample 
 * from is assumed to be a discrete probability distribution characterized 
 * by a collection of scenarios. Scenarios are assumed to be contained in a 
 * netCDF file and DiscreteScenarioSet gives method to deserialize the scenarios
 * from the netCDF file.
 * 
 * The deserialized scenario are stored in a boost::multi_array< double, 2 >.*/

class DiscreteScenarioSet : public ScenarioGenerator
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

    /// Container for the deserialized scenario pool
    /** Every scenario is assumed to have the same dimension. As the number of 
     * scenario becomes known whenever we deserialize the data, the scenario
     * pool is of known size at this point. Hence the choice to store it inside
     * a boost::multi_array.*/
    using DiscreteScenarioPool = boost::multi_array<double, 2>;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- PUBLIC VARIABLES ------------------------------*/
/*--------------------------------------------------------------------------*/
    /** @name Public variables
     *  @{ */

    /// Number of different scenarios in the scenario pool
    ScenarioIndex nbScenarios;

    /// Size of a scenario
    ScenarioSize scenarioSize;

    /// Container for Scenario-s
    DiscreteScenarioPool scenarioSet;

    // Pool size
    ScenarioIndex poolSize = 0;

    /// Indexes of the pool
    std::vector<ScenarioIndex> scenarioIndexes;

    /// Probabilities of scenarios
    std::vector<double> scenarioProbabilities;

    /// Random generator
    std::mt19937 rng;

    /// Current index in the pool
    ScenarioIndex currentScenarioIndex {0}; 

/** @} ---------------------------------------------------------------------*/
/*----------- CONSTRUCTING AND DESTRUCTING DiscreteScenarioSet -------------*/
/*--------------------------------------------------------------------------*/
    /** @name Constructing and destructing DiscreteScenarioSet
     *  @{ */

    DiscreteScenarioSet();

    /// deserialize a discrete distribution from a netCDF group
    /** Implementation of the "third-level" pure virtual function deserialize of 
     * ScenarioGenerator.h. Assumes that there is a two-dimensional variable
     *  \p Scenario contained inside a netCDF NcGroup. One dimension corresponds 
     * to the scenarioIndex and the other is the dimension component as a 
     * scenario is represented as a big vector in an euclidean space R^d. We 
     * deserialize it into a boost::multi_array< double, 2 > as the two 
     * dimensions become known once the file has been read.*/
    void deserialize(const netCDF::NcGroup &group) override;

    virtual ~DiscreteScenarioSet() = default;

    // Implementing pure virtual methods from ScenarioGenerator
    void set_seed(unsigned long seed) override;
    void init_random_pool(ScenarioIndex sampleSize) override;
    void init_representative_pool(ScenarioIndex sampleSize) override;
    Scenario get_current_scenario(void) override;
    double get_current_scenario_probability(void) override;
    bool next_scenario(void) override;
    ScenarioSize get_scenario_size(void) override;


/** @} ---------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

private:

/** @} ---------------------------------------------------------------------*/
/*-------------------- HELPER METHODS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
    /** @name helper methods of the class
     * Miscallenous functions
     * @{ */

    /// update the variable poolSize
    /** The first time a size for the pool has been given, update the 
     * variable poolSize accordingly. Also checks that the desired poolSize 
     * is less than the total number of available scenarios. */
   void update_poolSize(ScenarioIndex size) {
        if (size > nbScenarios)
        {
            throw std::out_of_range("The desired sample size is greater than "
            "the number of available number of different scenarios");
        }
        poolSize = size; 
   }

    /// Draw k elements among n
    /** The function generateRandomSubset draws k elements among n by use of 
     * the std::shuffle function and the internal rng. The chosen indexes are
     * moved into scenarioIndexes and these indexes characterize the 
     * scenarioPool used by the class DiscreteScenarioSet. */
    void generateRandomSubset(size_t n, size_t k)
    {
        if (k > n)
        {
            throw std::invalid_argument("k must be less or equal than n");
        }
        // elements = [1,2, ..., n]
        std::vector<ScenarioIndex> elements(n);
        std::iota(elements.begin(), elements.end(), 0);

        // Shuffle the elements randomly using our rng
        std::shuffle(elements.begin(), elements.end(), rng);

        // Move the first k elements into scenarioIndexes
        scenarioIndexes.resize(k);  
        std::move(elements.begin(), elements.end(), 
            scenarioIndexes.begin());
    }

/** @} ---------------------------------------------------------------------*/
/*------------------------- Kmeans clustering ------------------------------*/
/*--------------------------------------------------------------------------*/

    /** @name Kmeans clustering
     * 
     * Naive implementation from scratch of Lloyd's algorithm to solve the 
     * kmeans optimization problem. There are better existing librairies that 
     * are known to be robust and efficient to solve the kmeans optimization 
     * problem. The following implementation is simply here as a helper for
     * DiscreteScenarioSet to avoid additional dependencies with external
     * librairires.
     *
     * @{ */

    /// private type Point for Eigen::VectorXd
    /** For ease of linear algebra manipulations, we will work with 
     * Eigen::VectorXd, which stands morally for a scenario.
     * However, we will use the term Point to refer to a scenario represented 
     * as an Eigen::VectorXd.  
     * 
     * ?? Maybe a term like Eigen_scenario instead of Point is clearer ??
     * */
    using Point = Eigen::Map<Eigen::VectorXd>;

    /// Container type for the centers and the points to be clustered
    /** Lightweight container for a collection of scenarios represented as
     * Eigen::VectorXd. */
    using PoolMap = std::vector<Eigen::Map<Eigen::VectorXd>>;

    /// euclidean distance
    /** Computes the euclidean distance between two Eigen::VectorXd. */
    double euclideanDistance(const Eigen::VectorXd &vec1,
     const Eigen::VectorXd &vec2)
    {
        return (vec1 - vec2).norm();
    }

    /// (Costly) function to find the index of the nearest center to a point
    int nearestCenterIndex(const Eigen::VectorXd &point,
     const std::vector<Point> &centers)
    {
        double minDistance = std::numeric_limits<double>::max();
        int index = 0;
        for (size_t i = 0; i < centers.size(); i++)
        {
            double distance = euclideanDistance(point, centers[i]);
            if (distance < minDistance)
            {
                minDistance = distance;
                index = i;
            }
        }
        return index;
    }

    /// K-Means clustering function
    /** Naïve implementation of Lloyd's algorithm to solve the kmeans 
     * clustering problem. Return the centers of the clusters.
     * Could greatly benefit from better heuristics both in the centers 
     * initialization and in the stopping criteria.*/
    void kMeans(unsigned int k, PoolMap &pool) {

        size_t n = pool.size(); // n = nbScenarios
        std::vector<int> labels(n, 0);

        // Initialize centers randomly
        std::vector<Point> centers;
        for (int i = 0; i < k; i++)
        {
            centers.push_back(pool[rand() % n]);
        }

        bool changed;
        do
        {
            changed = false;

            // Assign points to the nearest center
            for (int i = 0; i < n; i++)
            {
                int newIndex = nearestCenterIndex(pool[i], centers);
                if (labels[i] != newIndex)
                {
                    labels[i] = newIndex;
                    changed = true;
                }
            }

            // Update centers by computing barycenter of each Voronoi cell
            std::vector<int> counts(k, 0);
            for (auto &center : centers)
            {
                center.setZero();
            }

            for (int i = 0; i < n; i++)
            {
                Point &center = centers[labels[i]];

                for (size_t j = 0; j < get_scenario_size(); j++)
                {
                    center[j] += pool[i][j];
                }
                counts[labels[i]]++;
            }

            for (int i = 0; i < k; i++)
            {
                for (size_t j = 0; j < get_scenario_size(); j++)
                {
                    centers[i][j] /= counts[i];
                }
            }
        } while (changed);

        // Replace the first size rows of scenarioSet with the cluster centers
        for (size_t i = 0; i < centers.size(); ++i) {
            std::copy(centers[i].data(), centers[i].data() +
                get_scenario_size(), scenarioSet[i].begin());
        }
    }

    // Macro for the factory
    SMSpp_insert_in_factory_h;

    /** @} */
};
} // end( namespace SMSpp_di_unipi_it )

#endif // __DiscreteScenarioSet

/*--------------------------------------------------------------------------*/
/*------------------- End file DiscreteScenarioSet.h -----------------------*/
/*--------------------------------------------------------------------------*/
