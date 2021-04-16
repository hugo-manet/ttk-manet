/// \ingroup vtk
/// \class ttkMultiTrackSource
/// \author Jonas Lukasczyk (jl@jluk.de)
/// \date 01.09.2019
///
/// This filter creates an MultiTrackSource with a specified radius, center, and
/// number of subdivisions.
///
/// \sa ttk::MultiTrackSource
/// \sa ttk::ttkAlgorithm

#pragma once

#include <random>
#include <vector>

// VTK Module
#include <ttkMultiTrackSourceModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkUnstructuredGrid.h>

// TTK Base Includes
#include <MultiTrackSource.h>

class vtkDataArray;

class TTKMULTITRACKSOURCE_EXPORT ttkMultiTrackSource
  : public ttkAlgorithm,
    protected ttk::MultiTrackSource {
private:
  int NumberOfLocations{3};
  double Radius{1};
  double StdDev{.25};
  int NbSmooth{3};
  double ZTranslation{.25};

  int TotalLength{20};
  int MinLength{3};
  int NbTrackPerLocation{2};
  double MassPortion{.5};

  std::default_random_engine generator;

  // method copied from ttkTrackingFromPD, but surely wrong
  vtkUnstructuredGrid *outputMesh_{nullptr};

public:
  static ttkMultiTrackSource *New();
  vtkTypeMacro(ttkMultiTrackSource, ttkAlgorithm);

  vtkSetMacro(NumberOfLocations, int);
  vtkGetMacro(NumberOfLocations, int);

  vtkSetMacro(NbTrackPerLocation, int);
  vtkGetMacro(NbTrackPerLocation, int);

  vtkSetMacro(Radius, double);
  vtkGetMacro(Radius, double);

  vtkSetMacro(StdDev, double);
  vtkGetMacro(StdDev, double);

  vtkSetMacro(NbSmooth, int);
  vtkGetMacro(NbSmooth, int);

  vtkSetMacro(ZTranslation, double);
  vtkGetMacro(ZTranslation, double);

  /*
  vtkSetMacro(, double);
  vtkGetMacro(, double);

  vtkSetMacro(, int);
  vtkGetMacro(, int);

  */

  vtkSetMacro(TotalLength, int);
  vtkGetMacro(TotalLength, int);

  vtkSetMacro(MinLength, int);
  vtkGetMacro(MinLength, int);

  vtkSetMacro(MassPortion, double);
  vtkGetMacro(MassPortion, double);

protected:
  ttkMultiTrackSource();
  ~ttkMultiTrackSource();

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  /* If number of timestep is even, peak value will be repeated,
   *  so that it is reached.
   * Coordinates std dev is with uniform distributions.
   * Coordinates smoothing generates nbSmooth more points,
   * and then takes the middle of segments nbSmooth times.
   */
  void addOneTrack(vtkUnstructuredGrid *dataset,
                   vtkPoints *points,
                   double center[3],
                   double baseStdDev[3],
                   int startTime,
                   int stopTimeExcl,
                   double peakPers,
                   int componentId,
                   vtkDoubleArray *costScalars,
                   vtkDoubleArray *persistenceScalars,
                   vtkDoubleArray *valueScalars,
                   vtkIntArray *matchingIdScalars,
                   vtkIntArray *lengthScalars,
                   vtkIntArray *timeScalars,
                   vtkIntArray *componentIds,
                   vtkIntArray *pointTypeScalars);

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
