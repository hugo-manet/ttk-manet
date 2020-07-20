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

      this->TimeLimit /= final_centroid.size();
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

    printMsg("Complete", 1, tm.getElapsedTime(), threadNumber_);
    return 1;
  }
} // namespace ttk
