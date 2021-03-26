#include "fuzzyDistanceMatrix.h"
#include <BottleneckDistance.h>
#include <TracksMatching.h>
#include <cstddef>
#include <limits>

namespace ttk {
  using namespace std;

  TracksMatching::TracksMatching() {
    // inherited from Debug: prefix will be printed at the beginning of every
    // msg
    this->setDebugMsgPrefix("TracksMatching");
  }

  double distPow(const TracksMatching::TimedPoint &ptA,
                 const TracksMatching::TimedPoint &ptB,
                 double p,
                 double tN,
                 double gL) {
    return sumPow(p, (ptA.val - ptB.val), (ptA.persistence - ptB.persistence),
                  gL * (ptA.x - ptB.x), gL * (ptA.y - ptB.y),
                  gL * (ptA.z - ptB.z), tN * (ptA.timeStep - ptB.timeStep));
  }

  Matrix getMatrixPow(const TracksMatching::Track &seqA,
                      const TracksMatching::Track &seqB,
                      double p,
                      double tN,
                      double gL) {
    size_t Ni = seqA.size(), Nj = seqB.size();
    Matrix ret(Ni, Nj);
    for(size_t x = 0; x < Ni; ++x)
      for(size_t y = 0; y < Nj; ++y)
        ret(x, y) = distPow(seqA[x], seqB[y], p, tN, gL);

    return ret;
  }

  vector<double> getCurvilinearDistPow(const TracksMatching::Track &seq,
                                       double p,
                                       double tN,
                                       double gL) {
    size_t N = seq.size();
    vector<double> ret(N - 1);
    for(size_t x = 0; x < N - 1; ++x)
      ret[x] = distPow(seq[x], seq[x + 1], p, tN, gL);

    return ret;
  }

  TracksMatching::TWED::TWED(Track *bidderTrack_,
                             Track *goodTrack_,
                             double p,
                             double lambda_p,
                             double tN,
                             double gL)
    : bidderTrack(bidderTrack_), goodTrack(goodTrack_) {
    if(bidderTrack->trackType != goodTrack->trackType) {
      this->value = std::numeric_limits<double>::infinity();
      return;
    }
    auto mat_p = getMatrixPow(*bidderTrack, *goodTrack, p, tN, gL);
    auto curvRow_p = getCurvilinearDistPow(*bidderTrack, p, tN, gL);
    auto curvCol_p = getCurvilinearDistPow(*goodTrack, p, tN, gL);
    double phi_p = pow(0.5, p - 1.);

    const size_t nRows = mat_p._nRows, nCols = mat_p._nCols;

    Matrix dyn(nRows, nCols, numeric_limits<double>::infinity());
    enum Direction { MATCH, DELROW, DELCOL, UNKNOWN };
    TMatrix<Direction> chosenDir(nRows, nCols, UNKNOWN);

    dyn(0, 0) = 0.;
    auto propagateTo = [&](size_t iRowStart, size_t jColStart, size_t deltaRow,
                           size_t deltaCol, Direction dir) {
      double newCost = dyn(iRowStart, jColStart);
      if(deltaCol == 0)
        newCost += lambda_p + phi_p * curvRow_p[iRowStart];
      else if(deltaRow == 0)
        newCost += lambda_p + phi_p * curvCol_p[jColStart];
      else
        newCost
          += mat_p(iRowStart + 1, jColStart + 1) + mat_p(iRowStart, jColStart);
      if(dyn(iRowStart + deltaRow, jColStart + deltaCol) > newCost) {
        dyn(iRowStart + deltaRow, jColStart + deltaCol) = newCost;
        chosenDir(iRowStart + deltaRow, jColStart + deltaCol) = dir;
      }
    };
    for(size_t iRow = 0; iRow < nRows - 1; ++iRow) {
      for(size_t jCol = 0; jCol < nCols - 1; ++jCol) {
        // Propagate in three directions
        propagateTo(iRow, jCol, 1, 0, DELROW);
        propagateTo(iRow, jCol, 0, 1, DELCOL);
        propagateTo(iRow, jCol, 1, 1, MATCH);
      }
      // Propagate along the last column
      propagateTo(iRow, nCols - 1, 1, 0, DELROW);
    }
    for(size_t jCol = 0; jCol < nCols - 1; ++jCol) {
      // Propagate along the last line
      propagateTo(nRows - 1, jCol, 0, 1, DELCOL);
    }

    this->value = pow(dyn(nRows - 1, nCols - 1), 1. / p);

    for(size_t iRow = nRows - 1, jCol = nCols - 1; iRow > 0 || jCol > 0;) {
      switch(chosenDir(iRow, jCol)) {
        case MATCH:
          matchedSegments.push_back({iRow - 1, jCol - 1});
          iRow--;
          jCol--;
          break;
        case DELROW:
          matchedSegments.push_back({iRow - 1, -1});
          iRow--;
          break;
        case DELCOL:
          matchedSegments.push_back({-1, jCol - 1});
          jCol--;
          break;
        case UNKNOWN:
          // unreachable
          break;
      };
    }
  }

  void TracksMatching::computeDistances(double p,
                                        double lambda_p,
                                        double timeNormalization,
                                        double geometricalLifting) {
    size_t nRows = bidders.size(), nCols = goods.size();
    distances = TMatrix<TWED>(nRows, nCols, TWED());
    for(size_t iRow = 0; iRow < nRows - 1; ++iRow)
      for(size_t jCol = 0; jCol < nCols - 1; ++jCol)
        distances(iRow, jCol) = TWED(&bidders[iRow], &goods[jCol], p, lambda_p,
                                     timeNormalization, geometricalLifting);
  }

  void TracksMatching::runMatching(double p) {
    /** TODO :
     *   - split by critical type
     *   - use algos that optimize for unbalanced matchings
     */
    size_t nBid = bidders.size(), nGood = goods.size();
    size_t totSize = nBid + nGood;
    vector<vector<double>> costMatrixForBundle(totSize);
    Munkres solver;

    for(size_t iRow = 0; iRow < nBid - 1; ++iRow) {
      for(size_t jCol = 0; jCol < nGood - 1; ++jCol)
        costMatrixForBundle[iRow].push_back(distances(iRow, jCol).value);
      double delCost = bidders[iRow].deletionCost(p);
      for(size_t jDel = 0; jDel < nBid - 1; ++jDel)
        costMatrixForBundle[iRow].push_back(delCost);
      for(size_t iDel = nBid; iDel < totSize - 1; ++iDel)
        costMatrixForBundle[iDel].push_back(delCost);
    }
    for(size_t iDel = nBid; iDel < totSize - 1; ++iDel)
      for(size_t jDel = nGood; jDel < totSize - 1; ++jDel)
        costMatrixForBundle[iDel].push_back(0.);

    vector<matchingTuple> solverOutput;
    solver.setInput(totSize, totSize, (void *)&costMatrixForBundle);
    solver.run(solverOutput);

    matchedTracks.clear();
    for(auto tup : solverOutput) {
      int iRow = get<0>(tup), jCol = get<1>(tup);
      if(iRow >= (int)nBid)
        iRow = -1;
      if(jCol >= (int)nGood)
        jCol = -1;
      if(iRow != -1 || jCol != -1)
        matchedTracks.push_back({iRow, jCol});
    }
  }
} // namespace ttk
