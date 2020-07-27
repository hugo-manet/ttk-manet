/// \ingroup base
/// \class ttk::ttkPersistenceDiagramBarycenter
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
/// \sa PersistenceTimeWarpClustering

#pragma once

#define BLocalMax ttk::CriticalType::Local_maximum
#define BLocalMin ttk::CriticalType::Local_minimum
#define BSaddle1 ttk::CriticalType::Saddle1
#define BSaddle2 ttk::CriticalType::Saddle2

// VTK includes -- to adapt
#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiBlockDataSetAlgorithm.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

// VTK Module
#include <ttkPersistenceTimeWarpClusteringModule.h>

// ttk code includes
#include <PersistenceDiagramBarycenter.h>
#include <PersistenceTimeWarpClustering.h>
//
//
#include <ttkAlgorithm.h>

// in this example, this wrapper takes a data-set on the input and produces a
// data-set on the output - to adapt.
// see the documentation of the vtkAlgorithm class to decide from which VTK
// class your wrapper should inherit.
class TTKPERSISTENCETIMEWARPCLUSTERING_EXPORT ttkPersistenceTimeWarpClustering
  : public ttkAlgorithm,
    protected ttk::PersistenceTimeWarpClustering {

public:
  // void setNumberOfInputsFromCommandLine(int number) {
  //   numberOfInputsFromCommandLine = number;
  //   SetNumberOfInputPorts(number);
  // }
  static ttkPersistenceTimeWarpClustering *New();

  vtkTypeMacro(ttkPersistenceTimeWarpClustering, ttkAlgorithm);

  // default ttk setters
  // void SetDebugLevel(int debugLevel) {
  //   setDebugLevel(debugLevel);
  //   Modified();
  //   needUpdate_ = true;
  // }

  // void SetThreads() {
  //   if(!UseAllCores)
  //     threadNumber_ = ThreadNumber;
  //   else {
  //     threadNumber_ = ttk::OsCall::getNumberOfCores();
  //   }
  //   Modified();
  //   needUpdate_ = true;
  // }

  /*void SetThreadNumber(int threadNumber){
    ThreadNumber = threadNumber;
    SetThreads();
  }*/

  // void SetUseAllCores(bool onOff) {
  //   UseAllCores = onOff;
  //   SetThreads();
  // }
  // // end of default ttk setters

  // set-getters macros to define from each variable you want to access from
  // the outside (in particular from paraview) - to adapt.

  // vtkSetMacro(ScalarField, std::string);
  // vtkGetMacro(ScalarField, std::string);

  vtkSetMacro(WassersteinMetric, int);
  vtkGetMacro(WassersteinMetric, int);

  void SetUseProgressive(int data) {
    UseProgressive = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(UseProgressive, int);

  void SetTimeLimit(double data) {
    TimeLimit = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(TimeLimit, double);

  // void SetThreadNumber(int data) {
  //   ThreadNumber = data;
  //   Modified();
  //   needUpdate_ = true;
  // }
  // vtkGetMacro(ThreadNumber, int);

  void SetAlpha(double data) {
    if(data > 0 && data <= 1) {

      Alpha = data;
    } else if(data > 1) {
      Alpha = 1;
    } else {
      Alpha = 0.001;
    }
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(Alpha, double);

  void SetAntiAlpha(double data) {
    double alpha = 1 - data;
    SetAlpha(alpha);
  }

  void SetDeltaLim(double data) {
    DeltaLim = data;
    Modified();
    needUpdate_ = true;
  }

  vtkGetMacro(DeltaLim, double);

  void SetLambda(double data) {
    Lambda = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(Lambda, double);

  void SetNumberOfClusters(int data) {
    NumberOfClusters = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(NumberOfClusters, int);

  void SetUseAccelerated(bool data) {
    UseAccelerated = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(UseAccelerated, bool);

  void SetUseKmeansppInit(bool data) {
    UseKmeansppInit = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(UseKmeansppInit, bool);

  void SetForceUseOfAlgorithm(bool data) {
    ForceUseOfAlgorithm = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(ForceUseOfAlgorithm, bool);

  void SetDeterministic(bool data) {
    Deterministic = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(Deterministic, bool);

  void SetPairTypeClustering(int data) {
    PairTypeClustering = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(PairTypeClustering, int);

  void SetSpacing(double spacing) {
    Spacing = spacing;
    oldSpacing = spacing;
    Modified();
  }
  vtkGetMacro(Spacing, double);

  void SetDisplayMethod(int displayMethod) {
    DisplayMethod = displayMethod;
    if(displayMethod == 0) { // compact display
      Spacing = 0;
    } else {
      Spacing = oldSpacing;
    }
    Modified();
  }

  vtkGetMacro(DisplayMethod, bool);

  void SetUseAdditionalPrecision(bool data) {
    UseAdditionalPrecision = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(UseAdditionalPrecision, bool);

  void SetDistanceWritingOptions(int data) {
    DistanceWritingOptions = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(DistanceWritingOptions, int);

  void SetUseInterruptible(bool data) {
    UseInterruptible = data;
    Modified();
    needUpdate_ = true;
  }
  vtkGetMacro(UseInterruptible, bool);

  void SetMethod(int method) {
    Method = method;
    needUpdate_ = true;
    Modified();
  }
  vtkGetMacro(Method, double);

  vtkSetMacro(NumberOfIterations, int);
  vtkGetMacro(NumberOfIterations, int);

  vtkSetMacro(DeletionCost, double);
  vtkGetMacro(DeletionCost, double);

protected:
  ttkPersistenceTimeWarpClustering();

  // using diagramType = ttk::DiagramTuple;

  using matchingType = std::tuple<ttk::SimplexId, ttk::SimplexId, double>;

  double getPersistenceDiagram(ttk::Diagram &diagram,
                               vtkUnstructuredGrid *CTPersistenceDiagram_);

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  vtkSmartPointer<vtkUnstructuredGrid> createMatchings();
  vtkSmartPointer<vtkMultiBlockDataSet> createOutputClusteredDiagrams();
  vtkSmartPointer<vtkMultiBlockDataSet> createOutputCentroids();
  vtkSmartPointer<vtkUnstructuredGrid> createOutputTimeWarp();
  vtkSmartPointer<vtkUnstructuredGrid> createOutputSlicePoints();

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  std::vector<ttk::DiagramCurve> intermediateDiagramsCurves_{};
  std::vector<std::vector<std::vector<matchingType>>> all_matchings_{};
  std::vector<std::vector<std::vector<TimeWarpTuple>>> time_warp_{};
  std::vector<ttk::DiagramCurve> final_centroid_{};
  std::vector<int> inv_clustering_{};

  // vtkUnstructuredGrid* output_clusters_;
  // vtkUnstructuredGrid* output_centroids_;

  double Spacing{1.0};
  int DisplayMethod{0};
  double oldSpacing{1.0};

  double max_dimension_total_{};
  int Method{0}; // 0 = progressive approach, 1 = Auction approach
  bool needUpdate_{true};

  // base code features
  // int doIt(const std::vector<vtkUnstructuredGrid *> &input,
  //          vtkUnstructuredGrid *outputClusters,
  //          vtkUnstructuredGrid *outputCentroids,
  //          vtkUnstructuredGrid *outputMatchings,
  //          int numInputs);

  // bool needsToAbort() override;

  // int updateProgress(const float &progress) override;
};
