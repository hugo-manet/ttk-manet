/// vTODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::DynamicTimeWarp
/// \author Hugo Manet <hugo.manet@ens.fr>
/// \date 2020-06-19
///
/// This module defines the %DynamicTimeWarp class that computes a discrete time warp
/// according to a given distance matrix
///
/// \b Related \b publication: \n TODO
/// 'DynamicTimeWarp'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2020.
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <Triangulation.h>

#include <boost/numeric/ublas/matrix.hpp>

namespace ttk {

  /**
   * The DynamicTimeWarp class provides methods to compute
   * an optimal warping path for a given distance matrix
   */
  class DynamicTimeWarp : virtual public Debug {

  public:
    DynamicTimeWarp() {
      this->setDebugMsgPrefix(
        "DynamicTimeWarp"); // inherited from Debug: prefix will be printed at the
      // beginning of every msg
    };
    ~DynamicTimeWarp(){};

    enum class Direction : std::uint8_t {
      DIR_SAME_ROW,
      DIR_BOTH,
      DIR_SAME_COL
    };

    /**
     * vTODO 3: Implmentation of the algorithm.
     *
     *         Note: If the algorithm requires a triangulation then this
     *               method must be called after the triangulation has been
     *               preconditioned for the upcoming operations.
     */

    std::vector<Direction> computeWarpingPath(
      const boost::numeric::ublas::matrix<double> &distanceMatrix,
      double DeletionCost) const;

  }; // DynamicTimeWarp class

} // namespace ttk
