///

#include <DynamicTimeWarp.h>
#include <limits>

using namespace std;
using namespace ttk;

vector<DynamicTimeWarp::Direction> DynamicTimeWarp::computeWarpingPath(
  const boost::numeric::ublas::matrix<double> &distanceMatrix,
  double DeletionCost) const {
  vector<DynamicTimeWarp::Direction> retVal;

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
  for(size_t iRow = 0; iRow < nRows - 1; ++iRow) {
    for(size_t jCol = 0; jCol < nCols - 1; ++jCol) {
      // Propagate in three directions
      if(dynCostPath(iRow + 1, jCol)
         < dynCostPath(iRow, jCol)
             + DeletionCost * distanceMatrix(iRow + 1, jCol)) {
        dynCostPath(iRow + 1, jCol)
          = dynCostPath(iRow, jCol)
            + DeletionCost * distanceMatrix(iRow + 1, jCol);
        pathDirection(iRow + 1, jCol)
          = DynamicTimeWarp::Direction::DIR_SAME_COL;
      }
      if(dynCostPath(iRow, jCol + 1)
         < dynCostPath(iRow, jCol)
             + DeletionCost * distanceMatrix(iRow, jCol + 1)) {
        dynCostPath(iRow, jCol + 1)
          = dynCostPath(iRow, jCol)
            + DeletionCost * distanceMatrix(iRow, jCol + 1);
        pathDirection(iRow, jCol + 1)
          = DynamicTimeWarp::Direction::DIR_SAME_ROW;
      }
      if(dynCostPath(iRow + 1, jCol + 1)
         < dynCostPath(iRow, jCol) + distanceMatrix(iRow + 1, jCol + 1)) {
        dynCostPath(iRow + 1, jCol + 1)
          = dynCostPath(iRow, jCol) + distanceMatrix(iRow + 1, jCol + 1);
        pathDirection(iRow + 1, jCol + 1)
          = DynamicTimeWarp::Direction::DIR_BOTH;
      }
    }
    // Propagate along the last column
    if(dynCostPath(iRow + 1, nCols - 1)
       < dynCostPath(iRow, nCols - 1)
           + DeletionCost * distanceMatrix(iRow + 1, nCols - 1)) {
      dynCostPath(iRow + 1, nCols - 1)
        = dynCostPath(iRow, nCols - 1)
          + DeletionCost * distanceMatrix(iRow + 1, nCols - 1);
      pathDirection(iRow + 1, nCols - 1)
        = DynamicTimeWarp::Direction::DIR_SAME_COL;
    }
  }
  for(size_t jCol = 0; jCol < nCols - 1; ++jCol) {
    // Propagate along the last line
    if(dynCostPath(nRows - 1, jCol + 1)
       < dynCostPath(nRows - 1, jCol)
           + DeletionCost * distanceMatrix(nRows - 1, jCol + 1)) {
      dynCostPath(nRows - 1, jCol + 1)
        = dynCostPath(nRows - 1, jCol)
          + DeletionCost * distanceMatrix(nRows - 1, jCol + 1);
      pathDirection(nRows - 1, jCol + 1)
        = DynamicTimeWarp::Direction::DIR_SAME_ROW;
    }
  }

  this->printMsg("Reconstructing path...");

  for(size_t iRow = nRows - 1, jCol = nCols - 1; iRow + jCol > 0;) {
    retVal.push_back(pathDirection(iRow, jCol));
    switch(pathDirection(iRow, jCol)) {
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
  }

  reverse(retVal.begin(), retVal.end());
  this->printMsg("Done reconstructing.");

  return retVal;
}
