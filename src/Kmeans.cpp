/*--------------------------------------------------------------------------*/
/*-------------------------- File Kmeans.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Concrete class kmeans that is a simple implementation from scratch of Lloyd's
 * algorithm to solve the kmeans optimization problem. There are better existing
 * librairies that are known to be robust and efficient to solve the kmeans
 * optimization problem. The following implementation is simply here as a 
 * default helper for DiscreteScenarioSet to avoid additional dependencies with
 * external librairires.
 *
 * In particular, to be coherent with the needs of DiscreteScenarioSet, the 
 * input of the class kmeans is expected to be a reference to a 
 * boost::multi_array< double, 2> (seen as a matrix) containing a finite number 
 * of Scenarios (vector of double of fixed size), one at each row. The output is
 * both the centers and the labels of each scenario which indicates their
 * respective center.
 *
 * \author Benoit Tran \n Dipartimento di Informatica \n Universita' di Pisa \n
 *
 * \copyright &copy; by Benoit Tran
 */

/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __Kmeans
#define __Kmeans
/* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SMSTypedefs.h"

#include <boost/multi_array.hpp>
#include <Eigen/Dense>
#include <iostream>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------ CLASS Kmeans ------------------------------*/
    class Kmeans
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

    /// Vectors of the euclidean space R^d are represented as Eigen vectors
    using Point = Eigen::Map<Eigen::VectorXd>;

    /// Container type for the centers and the points to be clustered
    using PoolMap = std::vector<Eigen::Map<Eigen::VectorXd>>;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- PUBLIC VARIABLES ------------------------------*/
/*--------------------------------------------------------------------------*/
    /** @name Public Types
     *  @{ */

    /// storing the dimension d of the scenarios in R^d
    size_t dimension;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- PUBLIC VARIABLES ------------------------------*/
/*--------------------------------------------------------------------------*/
        /** @name Public variables
         *  @{ */

        /// storing the input as Eigen::VectorXd
        PoolMap poolMap;

/** @} ---------------------------------------------------------------------*/
/*------------------ CONSTRUCTING AND DESTRUCTING kmeans -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing kmeans
 *  @{ */

/// constructing kmeans
/** Input should be a boost::multi_array< double, 2> which contains the 
 * scenarios. Initialization converts each row into an 
 * Eigen::VectorXd for easier linear algebra manipulations.
 * As boost::multi_array<double , 2> can be viewed as a matrix with row-major 
 * memory storage, if one knows the adress of the first element, then one knows 
 * the addresses of the whole row. We use Eigen::Map to avoid unecessary copies 
 * of the input. */
    Kmeans(boost::multi_array<double, 2> scenarioPool)
    {
        dimension = scenarioPool.shape()[1];
        // Use Eigen::Map to map each row to an Eigen::VectorXd
        for (size_t i = 0; i < scenarioPool.shape()[0]; i++)
        {
            poolMap.emplace_back(Eigen::Map<Eigen::VectorXd>(&scenarioPool[i][0], dimension));
        }
    }

    virtual ~Kmeans() = default;

    /// K-Means clustering function
    /** Naïve implementation of Lloyd's algorithm to solve the kmeans 
     * clustering problem. Return the centers of the clusters and labels of each point.
     * Could greatly benefit from better heuristics (in the centers initialization)
     * or in the stopping criteria.*/

    std::pair<std::vector<int>, std::vector<Point>> kMeans(int k, PoolMap &pool) {
        // std::srand(std::time(0));
        size_t n = pool.size();
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

                for (size_t j = 0; j < dimension; j++)
                {
                    center[j] += pool[i][j];
                }
                counts[labels[i]]++;
            }

            for (int i = 0; i < k; i++)
            {
                for (size_t j = 0; j < dimension; j++)
                {
                    centers[i][j] /= counts[i];
                }
            }
        } while (changed);

    return {labels, centers};    
}

/** @} ---------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

    private:

/** @} ---------------------------------------------------------------------*/
/*--------------------------- HELPER FUNCTIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

    /** @name helper functions
     * Helper functions for Kmeans.cpp
     * @{ */

    /// euclidean distance
    /**
     * Computes the euclidean distance between two Eigen::VectorXd
     */
    double euclideanDistance(const Eigen::VectorXd &vec1, const Eigen::VectorXd &vec2)
    {
        return (vec1 - vec2).norm();
    }

    /// (Costly) function to find the index of the nearest center to a point
    int nearestCenterIndex(const Eigen::VectorXd &point, const std::vector<Point> &centers)
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

    /** @} */
}; // class Kmeans

} // end( namespace SMSpp_di_unipi_it )

#endif // __Kmeans

/*--------------------------------------------------------------------------*/
/*------------------------- End file Kmeans.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
