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
#include <PersistenceDiagramBarycenter.h>
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
  class PersistenceTimeWarpClustering : virtual public Debug {

  public:
    PersistenceTimeWarpClustering() {
      this->setDebugMsgPrefix("PersistenceTimeWarpClustering");
    };

    ~PersistenceTimeWarpClustering(){};

    template <class dataType>
    int execute(
      std::vector<ttk::DiagramCurve> &intermediateDiagramCurves,
      ttk::DiagramCurve &final_centroid,
      std::vector<std::vector<std::vector<matchingTuple>>> &all_matchings);

    inline int setNumberOfInputs(int numberOfInputs) {
      numberOfInputs_ = numberOfInputs;
      // 			if(inputData_)
      // 			free(inputData_);
      // 			inputData_ = (void **) malloc(numberOfInputs*sizeof(void *));
      // 			for(int i=0 ; i<numberOfInputs ; i++){
      // 			inputData_[i] = NULL;
      // 			}
      return 0;
    }

    template <class dataType>
    static dataType abs(const dataType var) {
      return (var >= 0) ? var : -var;
    }

  protected:
    // Critical pairs used for clustering
    // 0:min-saddles ; 1:saddles-saddles ; 2:sad-max ; else : all

    int PairTypeClustering{-1};
    bool Deterministic{true};
    int WassersteinMetric{2};

    int numberOfInputs_{};
    bool UseProgressive{true};

    bool UseInterruptible{true};
    double Alpha{1.0};
    bool UseAdditionalPrecision{false};
    bool UseKmeansppInit{false};
    double DeltaLim{0.01};
    double Lambda{1.0};
    double TimeLimit{999999};

    int NumberOfClusters{1};
    bool UseAccelerated{false};

    int points_added_;
    int points_deleted_;

    // std::vector<BidderDiagram<void>> bidder_diagrams_;
    // std::vector<GoodDiagram<void>> barycenter_goods_;
  };

  template <class dataType>
  int PersistenceTimeWarpClustering::execute(
    std::vector<ttk::DiagramCurve> &intermediateDiagramCurves,
    ttk::DiagramCurve &final_centroid,
    std::vector<std::vector<std::vector<matchingTuple>>> &all_matchings) {

    Timer tm;
    {
      printMsg("Clustering " + std::to_string(numberOfInputs_) + " diagrams in "
               + std::to_string(NumberOfClusters) + " cluster(s).");
    }

    // No DTW, juste euclidian barycenter
    for(int iDiag = 0; iDiag < final_centroid.size(); ++iDiag) {
      std::vector<Diagram> slice;
      slice.reserve(numberOfInputs_);
      for(auto &curve : intermediateDiagramCurves) {

        slice.emplace_back(curve[iDiag]);
      }

      PersistenceDiagramBarycenter<double> persistenceDiagramsBarycenter;
      persistenceDiagramsBarycenter.setWasserstein("2");
      persistenceDiagramsBarycenter.setMethod(2);
      persistenceDiagramsBarycenter.setNumberOfInputs(numberOfInputs_);
      persistenceDiagramsBarycenter.setTimeLimit(TimeLimit);
      persistenceDiagramsBarycenter.setDeterministic(Deterministic);
      persistenceDiagramsBarycenter.setUseProgressive(UseProgressive);
      // persistenceDiagramsBarycenter.setDebugLevel(debugLevel_);
      // persistenceDiagramsBarycenter.setThreadNumber(threadNumber_);
      persistenceDiagramsBarycenter.setAlpha(Alpha);
      persistenceDiagramsBarycenter.setLambda(Lambda);
      // persistenceDiagramsBarycenter.setReinitPrices(ReinitPrices);
      // persistenceDiagramsBarycenter.setEpsilonDecreases(EpsilonDecreases);
      // persistenceDiagramsBarycenter.setEarlyStoppage(EarlyStoppage)

      persistenceDiagramsBarycenter.execute(
        slice, final_centroid[iDiag], all_matchings);
      ;
    }

    printMsg("Complete", 1, tm.getElapsedTime(), threadNumber_);
    return 1;
  }
} // namespace ttk
