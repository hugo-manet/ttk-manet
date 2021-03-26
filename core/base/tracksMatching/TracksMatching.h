/// TODO 1: Provide your information
///
/// \ingroup base
/// \class ttk::TracksMatching
/// \author Your Name Here <Your Email Address Here>
/// \date The Date Here.
///
/// This module defines the %TracksMatching class that computes...
///
/// \b Related \b publication: \n
/// 'TracksMatching'
/// Jonas Lukasczyk and Julien Tierny.
/// TTK Publications.
/// 2021.
///

#pragma once

// ttk common includes
#include "BottleneckDistance.h"
#include "DataTypes.h"
#include "fuzzyDistanceMatrix.h"
#include <Debug.h>
#include <Triangulation.h>
#include <limits>

namespace ttk {

  /**
   * The TracksMatching class provides methods to compute ...
   */
  class TracksMatching : virtual public Debug {
  public:
    using DataPoint = DecoratedDiagramTuple::DataPoint;
    struct TimedPoint : public DataPoint {
      TimedPoint(float timeStep_,
                 double val_,
                 float x_,
                 float y_,
                 float z_,
                 double persistence_)
        : DataPoint(val_, x_, y_, z_), timeStep(timeStep_),
          persistence(persistence_) {
      }
      float timeStep;
      double persistence;
    };
    struct Track : public std::vector<TimedPoint> {
      // = TWED_p cost of deletion: almost TWED_p to the same with pers := 0
      inline double deletionCost(double p) const {
        double ret = 0.;
        for(size_t i = 0; i < this->size(); ++i)
          ret += std::pow((*this)[i].persistence, p)
                 * (i != 0 && i != this->size() - 1 ? 2. : 1.);
        return std::pow(ret, 1. / p);
      }
      ttk::CriticalType trackType;
    };

    struct TWED {
      Track *bidderTrack;
      Track *goodTrack;
      double value;
      // (bidderId,goodId) ; -1 if deleted
      std::vector<std::pair<int, int>> matchedSegments;

      TWED()
        : bidderTrack(0), goodTrack(0),
          value(std::numeric_limits<double>::infinity()) {
      }
      TWED(Track *bidderTrack,
           Track *goodTrack,
           double p,
           double lambda_p,
           double timeNormalization,
           double geometricalLifting);

      inline int nbDeletion() const {
        return 2 * matchedSegments.size() - bidderTrack->size()
               - goodTrack->size();
      }
    };

    std::vector<Track> bidders, goods;
    TMatrix<TWED> distances;
    // (bidderId,goodId) ; -1 if deleted
    std::vector<std::pair<int, int>> matchedTracks;

  public:
    TracksMatching();
    void computeDistances(double p,
                          double lambda_p,
                          double timeNormalization,
                          double geometricalLifting);
    void runMatching(double p);
    inline std::vector<std::pair<int, int>> run(double p,
                                                double lambda_p,
                                                double timeNormalization,
                                                double geometricalLifting) {
      computeDistances(p, lambda_p, timeNormalization, geometricalLifting);
      runMatching(p);
      return matchedTracks;
    }
  }; // TracksMatching class

  template <typename... Args>
  double sumPow(double p, Args... vals) {
    return (std::pow(std::abs(vals), p) + ...);
  }
  /** This isn't the geom lifting of the distance in the (birth;death) space
   * (because we don't have access to it),
   * but almost : it's the (persistence;val) where val is either the birth or
   * the death. It's an equivalent norm though. Please don't complain.
   */
  double distPow(const TracksMatching::TimedPoint &ptA,
                 const TracksMatching::TimedPoint &ptB,
                 double p,
                 double timeNormalization,
                 double geometricalLifting);
  inline double dist(const TracksMatching::TimedPoint &ptA,
                     const TracksMatching::TimedPoint &ptB,
                     double p,
                     double timeNormalization,
                     double geometricalLifting) {
    return std::pow(
      distPow(ptA, ptB, p, timeNormalization, geometricalLifting), 1. / p);
  }
} // namespace ttk
