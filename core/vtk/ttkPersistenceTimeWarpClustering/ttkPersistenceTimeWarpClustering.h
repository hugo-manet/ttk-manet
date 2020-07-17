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
  static ttkPersistenceTimeWarpClustering *New();

  vtkTypeMacro(ttkPersistenceTimeWarpClustering, ttkAlgorithm);
  vtkSetMacro(WassersteinMetric, int);
  vtkGetMacro(WassersteinMetric, int);

  void SetUseProgressive(int data) {
    UseProgressive = data;
    Modified();
  }
  vtkGetMacro(UseProgressive, int);

  void SetTimeLimit(double data) {
    TimeLimit = data;
    Modified();
  }
  vtkGetMacro(TimeLimit, double);
  void SetAlpha(double data) {
    if(data > 0 && data <= 1) {

      Alpha = data;
    } else if(data > 1) {
      Alpha = 1;
    } else {
      Alpha = 0.001;
    }
    Modified();
  }
  vtkGetMacro(Alpha, double);

  void SetAntiAlpha(double data) {
    double alpha = 1 - data;
    SetAlpha(alpha);
  }

  void SetDeltaLim(double data) {
    DeltaLim = data;
    Modified();
  }

  vtkGetMacro(DeltaLim, double);

  void SetLambda(double data) {
    Lambda = data;
    Modified();
  }
  vtkGetMacro(Lambda, double);

  void SetNumberOfClusters(int data) {
    NumberOfClusters = data;
    Modified();
  }
  vtkGetMacro(NumberOfClusters, int);

  void SetUseAccelerated(bool data) {
    UseAccelerated = data;
    Modified();
  }
  vtkGetMacro(UseAccelerated, bool);

  void SetUseKmeansppInit(bool data) {
    UseKmeansppInit = data;
    Modified();
  }
  vtkGetMacro(UseKmeansppInit, bool);

  void SetDeterministic(bool data) {
    Deterministic = data;
    Modified();
  }
  vtkGetMacro(Deterministic, bool);

  void SetPairTypeClustering(int data) {
    PairTypeClustering = data;
    Modified();
  }
  vtkGetMacro(PairTypeClustering, int);

  void SetUseAdditionalPrecision(bool data) {
    UseAdditionalPrecision = data;
    Modified();
  }
  vtkGetMacro(UseAdditionalPrecision, bool);

  void SetUseInterruptible(bool data) {
    UseInterruptible = data;
    Modified();
  }
  vtkGetMacro(UseInterruptible, bool);

protected:
  ttkPersistenceTimeWarpClustering();

  using diagramType = ttk::DiagramTuple;

  using matchingType = std::tuple<ttk::SimplexId, ttk::SimplexId, double>;

  double getPersistenceDiagram(ttk::Diagram &diagram,
                               vtkUnstructuredGrid *CTPersistenceDiagram_);

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  vtkSmartPointer<vtkUnstructuredGrid> createMatchings();
  vtkSmartPointer<vtkUnstructuredGrid> createOutputClusteredDiagrams();
  vtkSmartPointer<vtkUnstructuredGrid> createOutputCentroids();

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

private:
  std::vector<ttk::DiagramCurve> intermediateDiagramsCurves_{};
  std::vector<std::vector<std::vector<matchingType>>> all_matchings_{};
  ttk::DiagramCurve final_centroid_{};
};
