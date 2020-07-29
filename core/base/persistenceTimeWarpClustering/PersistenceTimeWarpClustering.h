/// \ingroup base
/// \class ttk::PersistenceTimeWarpClustering
/// \author Jules Vidal <jules.vidal@lip6.fr>
/// \author Joseph Budin <joseph.budin@polytechnique.edu>
/// \date September 2019
///
/// \brief TTK processing package for the computation of Wasserstein barycenters
/// and K-Means clusterings of a set of persistence diagrams.
///
/// \b Related \b publication \n
/// "Progressive Wasserstein Barycenters of Persistence Diagrams" \n
/// Jules Vidal, Joseph Budin and Julien Tierny \n
/// Proc. of IEEE VIS 2019.\n
/// IEEE Transactions on Visualization and Computer Graphics, 2019.
///
/// \sa ttkPersistenceTimeWarpClustering

#pragma once

#ifndef BNodeType
#define BNodeType ttk::CriticalType
#define BLocalMax ttk::CriticalType::Local_maximum
#define BLocalMin ttk::CriticalType::Local_minimum
#define BSaddle1 ttk::CriticalType::Saddle1
#define BSaddle2 ttk::CriticalType::Saddle2
#define BIdVertex ttk::SimplexId
#endif

// base code includes
//
// #include <Wrapper.h>
//
// #include <PersistenceDiagram.h>
//
#include <DynamicTimeWarp.h>
#include <PersistenceDiagramClustering.h>
#include <PersistenceDiagramDistanceMatrix.h>
//
// #include <limits>
//
#include <PDClustering.h>
//

using namespace std;
using namespace ttk;

namespace ttk {

  using DiagramCurve = std::vector<ttk::Diagram>;
  using TimeWarpTuple = std::tuple<size_t, size_t, double>;
  class PersistenceTimeWarpClustering : public PersistenceDiagramClustering {
  protected:
    int NumberOfIterations{3};
    double DeletionCost{1.};
    int UseTWED{false};

  public:
    PersistenceTimeWarpClustering() {
      this->setDebugMsgPrefix("PersistenceTimeWarpClustering");
    };

    ~PersistenceTimeWarpClustering(){};

    template <class dataType>
    int executeTimeWarp(
      const std::vector<ttk::DiagramCurve> &intermediateDiagramCurves,
      ttk::DiagramCurve &final_centroid,
      std::vector<std::vector<std::vector<matchingTuple>>> &all_matchings,
      std::vector<std::vector<std::vector<TimeWarpTuple>>> &time_warp);
  };

  template <class dataType>
  int PersistenceTimeWarpClustering::executeTimeWarp(
    const std::vector<ttk::DiagramCurve> &intermediateDiagramCurves,
    ttk::DiagramCurve &final_centroid,
    std::vector<std::vector<std::vector<matchingTuple>>> &all_matchings,
    std::vector<std::vector<std::vector<TimeWarpTuple>>> &time_warp) {

    Timer tm;
#ifdef TTK_ENABLE_OPENMP
    const int svg_OpenMP_nested = omp_get_nested();
    omp_set_nested(1);
#endif // TTK_ENABLE_OPENMP

    const size_t nCurves = intermediateDiagramCurves.size();
    std::vector<std::vector<std::vector<std::pair<size_t, double>>>>
      matchedDiagrams(nCurves);
    for(auto &matchesForCurve : matchedDiagrams)
      matchesForCurve.assign(final_centroid.size(), {});
    {
      final_centroid = intermediateDiagramCurves[0];
      // list of all matched diagrams for centroid diagram
      for(int iIter = 0; iIter <= NumberOfIterations; ++iIter) {
        std::vector<std::vector<std::vector<std::pair<size_t, double>>>>
          oldMatchings(nCurves);
        for(auto &matchesForCurve : oldMatchings)
          matchesForCurve.assign(final_centroid.size(), {});
        std::swap(oldMatchings, matchedDiagrams);

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
        for(size_t jCurve = 0; jCurve < nCurves; ++jCurve) {
          Timer timerCurve;
          std::vector<Diagram> diagramColl;
          diagramColl.reserve(final_centroid.size()
                              + intermediateDiagramCurves[jCurve].size());
          diagramColl.insert(
            diagramColl.cend(), final_centroid.cbegin(), final_centroid.cend());
          diagramColl.insert(diagramColl.cend(),
                             intermediateDiagramCurves[jCurve].cbegin(),
                             intermediateDiagramCurves[jCurve].cend());
          const std::array<size_t, 2> sizes
            = {final_centroid.size(), intermediateDiagramCurves[jCurve].size()};

          printMsg("Computing distance matrix for centroid and curve "
                     + std::to_string(jCurve),
                   0, timerCurve.getElapsedTime(), threadNumber_);
          PersistenceDiagramDistanceMatrix distMatrixClass;
          distMatrixClass.setWasserstein(WassersteinMetric);
          distMatrixClass.setAlpha(Alpha);
          distMatrixClass.setLambda(Lambda);
          distMatrixClass.setDeltaLim(DeltaLim);
          // TODO ask Pierre if that's the best (without fuzzy matrices)
          // TODO put the enum in public visibility ^^
          // this is ConstraintType::RELATIVE_PERSISTENCE_PER_DIAG
          distMatrixClass.setConstraint(3);
          distMatrixClass.setThreadNumber(this->threadNumber_);

          auto distMatrix
            = distMatrixClass.execute(diagramColl, sizes, UseTWED);

          // TODO change dynTimeWarp input type to that less efficient one
          boost::numeric::ublas::matrix<double> realDistMatrix(
            sizes[0], sizes[1]);
          for(size_t i = 0; i < sizes[0]; ++i)
            for(size_t j = 0; j < sizes[1]; ++j)
              realDistMatrix(i, j) = distMatrix[i][j];

          printMsg("Time warping matrix for centroid and curve "
                     + std::to_string(jCurve),
                   0.5, timerCurve.getElapsedTime(), threadNumber_);
          std::vector<double> curvilinearDist;
          if(UseTWED) {
            curvilinearDist = std::move(distMatrix.back());
            distMatrix.pop_back();
          }
          // TODO parametrize from class param. We should also inherit DTW
          auto path = DynamicTimeWarp().computeWarpingPath(
            realDistMatrix, DeletionCost, UseTWED, curvilinearDist);

          double total_weight = realDistMatrix(0, 0);
          matchedDiagrams[jCurve][0].push_back({0, realDistMatrix(0, 0)});
          for(const auto &[dir, kDiagCentroid, lOther, w] : path) {
            total_weight += w;
            matchedDiagrams[jCurve][kDiagCentroid].push_back({lOther, w});
          }
          printMsg("Done with matrix for centroid and curve "
                     + std::to_string(jCurve) + ", distance from centroid "
                     + std::to_string(total_weight),
                   1, timerCurve.getElapsedTime(), threadNumber_);
        }
        size_t nbOfDifferentSlices = 0, nbOfDifferentMatch = 0;
        std::vector<bool> sliceChanged(final_centroid.size());
        for(size_t kDiag = 0; kDiag < final_centroid.size(); ++kDiag) {
          for(size_t jCurve = 0; jCurve < nCurves; ++jCurve) {
            const auto &oldie = oldMatchings[jCurve][kDiag];
            const auto &newbie = matchedDiagrams[jCurve][kDiag];
            size_t lOld = 0, lNew = 0;
            while(lOld < oldie.size() && lNew < newbie.size()) {
              if(oldie[lOld].first == newbie[lNew].first) {
                ++lOld;
                ++lNew;
              } else if(oldie[lOld].first < newbie[lNew].first) {
                sliceChanged[kDiag] = true;
                ++nbOfDifferentMatch;
                ++lOld;
              } else {
                sliceChanged[kDiag] = true;
                ++nbOfDifferentMatch;
                ++lNew;
              }
            }
            if(lOld < oldie.size() || lNew < newbie.size()) {
              sliceChanged[kDiag] = true;
              nbOfDifferentMatch += oldie.size() - lOld + newbie.size() - lNew;
            }
          }
          if(sliceChanged[kDiag])
            ++nbOfDifferentSlices;
        }
        if(iIter == NumberOfIterations || nbOfDifferentSlices == 0)
          break;
        const int svg_DebugLevel = this->debugLevel_;
        this->setDebugLevel(1);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
        for(size_t kDiag = 0; kDiag < final_centroid.size(); ++kDiag) {
          if(!sliceChanged[kDiag])
            continue;
          std::vector<Diagram> slice;
          for(size_t jCurve = 0; jCurve < nCurves; ++jCurve) {
            for(auto [lOther, w] : matchedDiagrams[jCurve][kDiag])
              slice.emplace_back(intermediateDiagramCurves[jCurve][lOther]);
          }

          std::vector<std::vector<std::vector<matchingTuple>>> temp_matchings;
          std::vector<Diagram> solo_centroid(1);
          this->execute<dataType>(slice, solo_centroid, temp_matchings);
          final_centroid[kDiag] = std::move(solo_centroid[0]);
          solo_centroid.clear();
/* No need to copy for now, it's just an overhead
          if(iIter >= NumberOfIterations - 1) {
            all_matchings.reserve(all_matchings.size() + temp_matchings.size());
            std::move(std::begin(temp_matchings), std::end(temp_matchings),
                      std::back_inserter(all_matchings));
            temp_matchings.clear();
          } // */
        }
        this->setDebugLevel(svg_DebugLevel);
        printMsg("Completed iteration with "
                   + std::to_string(nbOfDifferentMatch) + " changes in "
                   + std::to_string(nbOfDifferentSlices) + " slices",
                 iIter / (double)NumberOfIterations, tm.getElapsedTime(),
                 threadNumber_);
      }
    }
#ifdef TTK_ENABLE_OPENMP
    omp_set_nested(svg_OpenMP_nested);
#endif // TTK_ENABLE_OPENMP

    printMsg("Completed all iterations", 1, tm.getElapsedTime(), threadNumber_);
    time_warp.emplace_back(nCurves); // Only one cluster, all curves in
    for(size_t jCurve = 0; jCurve < nCurves; ++jCurve) {
      for(size_t kDiag = 0; kDiag < final_centroid.size(); ++kDiag) {
        for(auto [lOther, w] : matchedDiagrams[jCurve][kDiag]) {
          time_warp[0][jCurve].emplace_back(kDiag, lOther, w);
        }
      }
    }
    return 1;
  }
} // namespace ttk
