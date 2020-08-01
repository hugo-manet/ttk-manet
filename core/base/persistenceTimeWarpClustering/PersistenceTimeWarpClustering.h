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

template <class A, class B>
std::ostream &operator<<(std::ostream &os, const std::tuple<A, A, B> &p) {
  return os << '{' << std::get<0>(p) << ',' << std::get<1>(p) << ','
            << std::get<2>(p) << '}';
}

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
    std::vector<size_t> offsetForCurve = {0};
    std::vector<size_t> nDiagOfCurve(nCurves);
    for(size_t jCurve = 0; jCurve < nCurves; ++jCurve) {
      nDiagOfCurve[jCurve] = intermediateDiagramCurves[jCurve].size();
      offsetForCurve.push_back(offsetForCurve.back() + nDiagOfCurve[jCurve]);
    }
    std::vector<std::vector<std::tuple<int, int, double>>> matchGraph(
      offsetForCurve.back() + final_centroid.size());
    std::vector<std::vector<std::vector<std::pair<int, double>>>>
      matchedDiagrams(nCurves);
    std::vector<double> totalEnergy;
    for(auto &matchesForCurve : matchedDiagrams)
      matchesForCurve.assign(final_centroid.size(), {});
    {
      final_centroid = intermediateDiagramCurves[0];
      std::vector<bool> sliceChanged(final_centroid.size(), true);
      std::vector<boost::numeric::ublas::matrix<double>> realDistMatrix(
        nCurves);
      for(size_t jCurve = 0; jCurve < nCurves; ++jCurve)
        realDistMatrix[jCurve].resize(
          final_centroid.size(), intermediateDiagramCurves[jCurve].size());

      for(int iIter = 0; iIter <= NumberOfIterations; ++iIter) {
        std::vector<std::vector<std::tuple<int, int, double>>> oldMatchGraph(
          matchGraph.size());
        std::swap(oldMatchGraph, matchGraph);
        for(auto &matchesForCurve : matchedDiagrams)
          matchesForCurve.assign(final_centroid.size(), {});
        totalEnergy.push_back(0);

#ifdef TTK_ENABLE_OPENMP
//#pragma omp parallel for schedule(dynamic) num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
        for(size_t jCurve = 0; jCurve < nCurves; ++jCurve) {
          Timer timerCurve;
          std::vector<Diagram> diagramColl;
          diagramColl.reserve(final_centroid.size()
                              + intermediateDiagramCurves[jCurve].size());
          for(size_t kDiag = 0; kDiag < final_centroid.size(); ++kDiag)
            if(sliceChanged[kDiag])
              diagramColl.push_back(final_centroid[kDiag]);
          const std::array<size_t, 2> sizes
            = {diagramColl.size(), intermediateDiagramCurves[jCurve].size()};
          diagramColl.insert(diagramColl.cend(),
                             intermediateDiagramCurves[jCurve].cbegin(),
                             intermediateDiagramCurves[jCurve].cend());

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

          for(size_t kDiag = 0, kComputed = 0; kDiag < sizes[0]; ++kDiag) {
            if(sliceChanged[kDiag]) {
              for(size_t lOther = 0; lOther < sizes[1]; ++lOther)
                realDistMatrix[jCurve](kDiag, lOther)
                  = distMatrix[kComputed][lOther];
              ++kComputed;
            }
          }

          printMsg("Time warping matrix for centroid and curve "
                     + std::to_string(jCurve),
                   0.5, timerCurve.getElapsedTime(), threadNumber_);
          std::vector<double> curvilinearDist;
          if(UseTWED) {
            curvilinearDist = std::move(distMatrix.back());
            distMatrix.pop_back();
          }
          DynamicTimeWarp dtwClass;
          dtwClass.UseTWED = UseTWED;
          dtwClass.DeletionCost = DeletionCost;
          dtwClass.MetricExponent = 2.;
          // TODO parametrize from class param. We should also inherit DTW
          auto path = dtwClass.computeWarpingPath(
            realDistMatrix[jCurve], curvilinearDist);

          double total_weight = realDistMatrix[jCurve](0, 0);
          if(!UseTWED) {
            matchedDiagrams[jCurve][0].push_back(
              {0, realDistMatrix[jCurve](0, 0)});
            matchGraph[offsetForCurve.back()].push_back(
              {jCurve, 0, total_weight});
            matchGraph[offsetForCurve[jCurve]].push_back({-1, 0, total_weight});
          }
          for(const auto &[dir, kDiagCentroid, lOther, w] : path) {
            total_weight += w;
            if(UseTWED) {
              switch(dir) {
                case DynamicTimeWarp::Direction::DIR_SAME_COL:
                  matchedDiagrams[jCurve][kDiagCentroid].push_back({-1, w});
                  matchGraph[offsetForCurve.back() + kDiagCentroid - 1]
                    .push_back({-1, kDiagCentroid, w});
                  matchGraph[offsetForCurve.back() + kDiagCentroid].push_back(
                    {-1, kDiagCentroid - 1, w});
                  break;
                case DynamicTimeWarp::Direction::DIR_BOTH:
                  matchedDiagrams[jCurve][kDiagCentroid].push_back({lOther, w});
                  matchGraph[offsetForCurve[jCurve] + lOther].push_back(
                    {-1, kDiagCentroid, w});
                  matchGraph[offsetForCurve.back() + kDiagCentroid].push_back(
                    {jCurve, lOther, w});
                  matchGraph[offsetForCurve[jCurve] + lOther - 1].push_back(
                    {-1, kDiagCentroid - 1, w});
                  matchGraph[offsetForCurve.back() + kDiagCentroid - 1]
                    .push_back({jCurve, lOther - 1, w});
                  break;
                case DynamicTimeWarp::Direction::DIR_SAME_ROW:
                  matchGraph[offsetForCurve[jCurve] + lOther].push_back(
                    {jCurve, lOther - 1, w});
                  matchGraph[offsetForCurve[jCurve] + lOther - 1].push_back(
                    {jCurve, lOther, w});
                  break;
              }
            } else {
              matchedDiagrams[jCurve][kDiagCentroid].push_back({lOther, w});
              matchGraph[offsetForCurve[jCurve] + lOther].push_back(
                {-1, kDiagCentroid, w});
              matchGraph[offsetForCurve.back() + kDiagCentroid].push_back(
                {jCurve, lOther, w});
            }
          }
          printMsg("Done with matrix for centroid and curve "
                     + std::to_string(jCurve) + ", distance from centroid "
                     + std::to_string(total_weight),
                   1, timerCurve.getElapsedTime(), threadNumber_);
          totalEnergy.back() += total_weight;
        }
        printMsg("Total energy for step " + std::to_string(iIter) + ": "
                   + std::to_string(totalEnergy.back()),
                 1, tm.getElapsedTime(), threadNumber_);
        size_t nbOfDifferentSlices = 0, nbOfDifferentMatch = 0;
        sliceChanged.assign(final_centroid.size(), false);
        for(size_t kDiag = 0; kDiag < final_centroid.size(); ++kDiag) {
          const auto &oldie = oldMatchGraph[offsetForCurve.back() + kDiag];
          const auto &newbie = matchGraph[offsetForCurve.back() + kDiag];
          size_t lOld = 0, lNew = 0;
          while(lOld < oldie.size() && lNew < newbie.size()) {
            if(std::get<0>(oldie[lOld]) == std::get<0>(newbie[lNew])
               && std::get<1>(oldie[lOld]) == std::get<1>(newbie[lNew])) {
              if(std::get<2>(oldie[lOld]) != std::get<2>(newbie[lNew]))
                if(this->debugLevel_ > 3)
                  std::cout << oldie[lOld] << " != " << newbie[lNew]
                            << " in slice " << kDiag << std::endl;
              ++lOld;
              ++lNew;
            } else if(std::get<0>(oldie[lOld]) < std::get<0>(newbie[lNew])
                      || (std::get<0>(oldie[lOld]) == std::get<0>(newbie[lNew])
                          && std::get<1>(oldie[lOld])
                               < std::get<1>(newbie[lNew]))) {
              if(this->debugLevel_ > 3)
                std::cout << "Removed " << oldie[lOld] << " from slice "
                          << kDiag << std::endl;
              sliceChanged[kDiag] = true;
              ++nbOfDifferentMatch;
              ++lOld;
            } else {
              if(this->debugLevel_ > 3)
                std::cout << "Added " << newbie[lNew] << " to slice " << kDiag
                          << std::endl;
              sliceChanged[kDiag] = true;
              ++nbOfDifferentMatch;
              ++lNew;
            }
          }
          while(lOld < oldie.size()) {
            if(this->debugLevel_ > 3)
              std::cout << "Removed " << oldie[lOld] << " from slice " << kDiag
                        << std::endl;
            sliceChanged[kDiag] = true;
            ++nbOfDifferentMatch;
            ++lOld;
          }
          while(lNew < newbie.size()) {
            if(this->debugLevel_ > 3)
              std::cout << "Added " << newbie[lNew] << " to slice " << kDiag
                        << std::endl;
            sliceChanged[kDiag] = true;
            ++nbOfDifferentMatch;
            ++lNew;
          }

          if(sliceChanged[kDiag]) {
            ++nbOfDifferentSlices;
            double oldWeightOfSlice = 0., newWeightOfSlice = 0.;
            for(auto [jC, lO, oW] : oldMatchGraph[kDiag]) {
              oldWeightOfSlice += oW;
              if(this->debugLevel_ > 3)
                std::cout << "old (" << jC << ',' << lO << "," << oW << ")"
                          << std::endl;
            }
            for(auto [jC, lO, nW] : matchGraph[kDiag]) {
              newWeightOfSlice += nW;
              if(this->debugLevel_ > 3)
                std::cout << "new (" << jC << ',' << lO << "," << nW << ")"
                          << std::endl;
            }
            std::cout << "Slice " << kDiag << " changed from "
                      << oldWeightOfSlice << " to " << newWeightOfSlice
                      << std::endl;
          }
        }
        if(iIter == NumberOfIterations || nbOfDifferentSlices == 0)
          break;
        const int svg_DebugLevel = this->debugLevel_;
        this->setDebugLevel(1);
        auto recomputeDiagram = [&](size_t kDiag) {
          if(!sliceChanged[kDiag])
            return;
          std::vector<Diagram> slice;
          for(size_t jCurve = 0; jCurve < nCurves; ++jCurve) {
            if(UseTWED) {
              for(auto [lOther, w] : matchedDiagrams[jCurve][kDiag]) {
                if(svg_DebugLevel > 3)
                  std::cout << "Found " << lOther << " for " << jCurve << ","
                            << kDiag << " in loop A" << std::endl;
                if(lOther != -1) {
                  // Need to emplace twice to optimize the real distance
                  slice.emplace_back(intermediateDiagramCurves[jCurve][lOther]);
                  slice.emplace_back(intermediateDiagramCurves[jCurve][lOther]);
                } else
                  slice.emplace_back(final_centroid[kDiag - 1]);
              }
              if(kDiag + 1 < final_centroid.size()) {
                for(auto [lOther, w] : matchedDiagrams[jCurve][kDiag + 1]) {
                  if(svg_DebugLevel > 3)
                    std::cout << "Found " << lOther << " for " << jCurve << ","
                              << kDiag << " in loop B" << std::endl;
                  if(lOther != -1) {
                    // They matched on the following step, so we're also matched
                    slice.emplace_back(
                      intermediateDiagramCurves[jCurve][lOther - 1]);
                    slice.emplace_back(
                      intermediateDiagramCurves[jCurve][lOther - 1]);
                  } else
                    slice.emplace_back(final_centroid[kDiag + 1]);
                }
              }
            } else
              for(auto [lOther, w] : matchedDiagrams[jCurve][kDiag])
                slice.emplace_back(intermediateDiagramCurves[jCurve][lOther]);
          }

          if(svg_DebugLevel > 3)
            std::cout << "Slice is now of length " << slice.size() << std::endl;
          std::vector<std::vector<std::vector<matchingTuple>>> temp_matchings;
          std::vector<Diagram> solo_centroid(1);
          this->execute<dataType>(slice, solo_centroid, temp_matchings);
          final_centroid[kDiag] = std::move(solo_centroid[0]);
          solo_centroid.clear();
          /* No need to copy for now, it's just an overhead
                    if(iIter >= NumberOfIterations - 1) {
                      all_matchings.reserve(all_matchings.size() +
             temp_matchings.size()); std::move(std::begin(temp_matchings),
             std::end(temp_matchings), std::back_inserter(all_matchings));
                      temp_matchings.clear();
                    } // */
        };
        // TODO maybe lock the dependencies with mutexes to mix parallels ?
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
        for(size_t kDiag = 0; kDiag < final_centroid.size(); kDiag += 2)
          recomputeDiagram(kDiag);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif // TTK_ENABLE_OPENMP
        for(size_t kDiag = 1; kDiag < final_centroid.size(); kDiag += 2)
          recomputeDiagram(kDiag);
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
    std::cout << "Energy sequence : {";
    for(size_t i = 0; i < totalEnergy.size(); ++i)
      std::cout << ", " << i << ": " << totalEnergy[i];
    std::cout << "}" << std::endl;
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
