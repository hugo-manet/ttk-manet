///

#include <DynamicTimeWarp.h>
#include <limits>

using namespace std;
using namespace ttk;

vector<tuple<DynamicTimeWarp::Direction, size_t, size_t, double>>
  DynamicTimeWarp::computeWarpingPath(
    const boost::numeric::ublas::matrix<double> &distanceMatrix,
    double DeletionCost) const {
  vector<tuple<DynamicTimeWarp::Direction, size_t, size_t, double>> retVal;

  size_t nRows = distanceMatrix.size1(), nCols = distanceMatrix.size2();
  boost::numeric::ublas::matrix<double> dynCostPath(
    nRows, nCols, numeric_limits<double>::infinity());
  boost::numeric::ublas::matrix<DynamicTimeWarp::Direction> pathDirection(
    nRows, nCols, DynamicTimeWarp::Direction::DIR_BOTH);

  this->printMsg("Solving the dynamic algorithm by propagation...");
  /** We solve the dynamic problem by propagation
   *  We locally update the neighbourgs if we have a better path
   *  And if so, we update the preferred direction
   *  The propagation on the last column and last line are different to avoid
   * going out of the table
   */
  dynCostPath(0, 0) = distanceMatrix(0, 0);
  auto propagateTo
    = [&](size_t iRowStart, size_t jColStart, size_t deltaRow, size_t deltaCol,
          DynamicTimeWarp::Direction dir, double multiplier) {
        if(dynCostPath(iRowStart + deltaRow, jColStart + deltaCol)
           > dynCostPath(iRowStart, jColStart)
               + multiplier
                   * (distanceMatrix(iRowStart + deltaRow, jColStart + deltaCol)
                      + distanceMatrix(iRowStart, jColStart))) {
          dynCostPath(iRowStart + deltaRow, jColStart + deltaCol)
            = dynCostPath(iRowStart, jColStart)
              + multiplier
                  * (distanceMatrix(iRowStart + deltaRow, jColStart + deltaCol)
                     + distanceMatrix(iRowStart, jColStart));
          pathDirection(iRowStart + deltaRow, jColStart + deltaCol) = dir;
        }
      };
  for(size_t iRow = 0; iRow < nRows - 1; ++iRow) {
    for(size_t jCol = 0; jCol < nCols - 1; ++jCol) {
      // Propagate in three directions
      propagateTo(iRow, jCol, 1, 0, DynamicTimeWarp::Direction::DIR_SAME_COL,
                  DeletionCost);
      propagateTo(iRow, jCol, 0, 1, DynamicTimeWarp::Direction::DIR_SAME_ROW,
                  DeletionCost);
      propagateTo(iRow, jCol, 1, 1, DynamicTimeWarp::Direction::DIR_BOTH, 1);
    }
    // Propagate along the last column
    propagateTo(iRow, nCols - 1, 1, 0, DynamicTimeWarp::Direction::DIR_SAME_COL,
                DeletionCost);
  }
  for(size_t jCol = 0; jCol < nCols - 1; ++jCol) {
    // Propagate along the last line
    propagateTo(nRows - 1, jCol, 0, 1, DynamicTimeWarp::Direction::DIR_SAME_ROW,
                DeletionCost);
  }

  this->printMsg("Reconstructing path...");

  for(size_t iRow = nRows - 1, jCol = nCols - 1; iRow + jCol > 0;) {
    auto dir = pathDirection(iRow, jCol);
    double multiplier
      = (dir == DynamicTimeWarp::Direction::DIR_BOTH) ? 1 : DeletionCost;
    // pushed with a wrong (incomplete) weight, completed after
    retVal.push_back(
      {dir, iRow, jCol, distanceMatrix(iRow, jCol) * multiplier});
    switch(dir) {
      case DynamicTimeWarp::Direction::DIR_BOTH:
        --iRow;
        --jCol;
        break;
      case DynamicTimeWarp::Direction::DIR_SAME_ROW:
        --jCol;
        break;
      case DynamicTimeWarp::Direction::DIR_SAME_COL:
        --iRow;
        break;
    }
    // complete
    get<3>(*(retVal.rbegin())) += distanceMatrix(iRow, jCol) * multiplier;
  }

  reverse(retVal.begin(), retVal.end());
  this->printMsg("Done reconstructing.");

  return retVal;
}
