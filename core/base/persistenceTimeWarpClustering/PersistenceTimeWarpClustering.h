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
  class PersistenceTimeWarpClustering : public PersistenceDiagramClustering {

  public:
    PersistenceTimeWarpClustering() {
      this->setDebugMsgPrefix("PersistenceTimeWarpClustering");
    };

    ~PersistenceTimeWarpClustering(){};

    template <class dataType>
    int executeTimeWarp(
      std::vector<ttk::DiagramCurve> &intermediateDiagramCurves,
      ttk::DiagramCurve &final_centroid,
      std::vector<std::vector<std::vector<matchingTuple>>> &all_matchings);
  };

  template <class dataType>
  int PersistenceTimeWarpClustering::executeTimeWarp(
    std::vector<ttk::DiagramCurve> &intermediateDiagramCurves,
    ttk::DiagramCurve &final_centroid,
    std::vector<std::vector<std::vector<matchingTuple>>> &all_matchings) {

    Timer tm;

    const size_t nCurves = intermediateDiagramCurves.size();

    if(false) {
      this->TimeLimit /= final_centroid.size();
      // No DTW, juste euclidian barycenter
      for(int iDiag = 0; iDiag < final_centroid.size(); ++iDiag) {
        std::vector<Diagram> slice;
        slice.reserve(numberOfInputs_);
        for(auto &curve : intermediateDiagramCurves) {
          slice.emplace_back(curve[iDiag]);
        }
        numberOfInputs_ = slice.size();
        {
          printMsg("Clustering " + std::to_string(numberOfInputs_)
                   + " diagrams in " + std::to_string(NumberOfClusters)
                   + " cluster(s).");
        }

        std::vector<std::vector<std::vector<matchingTuple>>> temp_matchings;
        std::vector<Diagram> solo_centroid(1);
        this->execute<dataType>(slice, solo_centroid, temp_matchings);
        final_centroid[iDiag] = std::move(solo_centroid[0]);
        solo_centroid.clear();
        all_matchings.reserve(all_matchings.size() + temp_matchings.size());
        std::move(std::begin(temp_matchings), std::end(temp_matchings),
                  std::back_inserter(all_matchings));
        temp_matchings.clear();
      }
    } else {
      final_centroid = intermediateDiagramCurves[0];
      // TODO param
      size_t nbIterMax = 3;
      // list of all matched diagrams for centroid diagram
      for(size_t iIter = 0; iIter < nbIterMax; ++iIter) {
        std::vector<std::vector<std::vector<size_t>>> matchedDiagrams(
          final_centroid.size());
        for(auto &matchesOfDiag : matchedDiagrams)
          matchesOfDiag.assign(nCurves, {});

        for(size_t jCurve = 0; jCurve < nCurves; ++jCurve) {
          std::vector<Diagram> diagramColl = final_centroid;
          diagramColl.insert(diagramColl.end(),
                             intermediateDiagramCurves[jCurve].begin(),
                             intermediateDiagramCurves[jCurve].end());
          const std::array<size_t, 2> sizes
            = {final_centroid.size(), intermediateDiagramCurves[jCurve].size()};

          PersistenceDiagramDistanceMatrix distMatrixClass;
          distMatrixClass.setWasserstein(WassersteinMetric);
          distMatrixClass.setAlpha(Alpha);
          distMatrixClass.setLambda(Lambda);
          distMatrixClass.setDeltaLim(DeltaLim);
          // TODO ask Pierre if that's the best (without fuzzy matrices)
          // TODO put the enum in public visibility ^^
          // this is ConstraintType::RELATIVE_PERSISTENCE_PER_DIAG
          distMatrixClass.setConstraint(3);

          auto distMatrix = distMatrixClass.execute(diagramColl, sizes);

          // TODO change dynTimeWarp input type to that less efficient one
          boost::numeric::ublas::matrix<double> realDistMatrix(
            sizes[0], sizes[1]);
          for(size_t i = 0; i < sizes[0]; ++i)
            for(size_t j = 0; j < sizes[1]; ++j)
              realDistMatrix(i, j) = distMatrix[i][j];

          // TODO parametrize from class param. We should also inherit DTW
          auto path
            = DynamicTimeWarp().computeWarpingPath(realDistMatrix, 1, false);

          matchedDiagrams[0][jCurve].push_back(0);
          for(const auto &[dir, iCentroid, kOther, w] : path) {
            matchedDiagrams[iCentroid][jCurve].push_back(kOther);
          }
        }
        for(int iDiag = 0; iDiag < final_centroid.size(); ++iDiag) {
          std::vector<Diagram> slice;
          for(size_t jCurve = 0; jCurve < nCurves; ++jCurve) {
            for(auto &kOther : matchedDiagrams[iDiag][jCurve])
              slice.emplace_back(intermediateDiagramCurves[jCurve][kOther]);
          }
          numberOfInputs_ = slice.size();
          {
            printMsg("Clustering " + std::to_string(numberOfInputs_)
                     + " diagrams in " + std::to_string(NumberOfClusters)
                     + " cluster(s).");
          }

          std::vector<std::vector<std::vector<matchingTuple>>> temp_matchings;
          std::vector<Diagram> solo_centroid(1);
          this->execute<dataType>(slice, solo_centroid, temp_matchings);
          final_centroid[iDiag] = std::move(solo_centroid[0]);
          solo_centroid.clear();

          if(iIter >= nbIterMax - 1) {
            all_matchings.reserve(all_matchings.size() + temp_matchings.size());
            std::move(std::begin(temp_matchings), std::end(temp_matchings),
                      std::back_inserter(all_matchings));
            temp_matchings.clear();
          }
        }
      }
    }

    printMsg("Complete", 1, tm.getElapsedTime(), threadNumber_);
    return 1;
  }
} // namespace ttk
