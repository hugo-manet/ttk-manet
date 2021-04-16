#include <algorithm>
#include <random>
#include <set>
#include <ttkMultiTrackSource.h>

#include <ttkUtils.h>

#include <vtkInformation.h>
#include <vtkObjectFactory.h>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

vtkStandardNewMacro(ttkMultiTrackSource);

ttkMultiTrackSource::ttkMultiTrackSource() {
  this->SetNumberOfInputPorts(0);
  this->SetNumberOfOutputPorts(1);
  generator.seed(time(NULL));
}
ttkMultiTrackSource::~ttkMultiTrackSource() {
  if(outputMesh_)
    outputMesh_->Delete();
}

int ttkMultiTrackSource::FillInputPortInformation(int port,
                                                  vtkInformation *info) {
  return 0;
}

int ttkMultiTrackSource::FillOutputPortInformation(int port,
                                                   vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  else
    return 0;
  return 1;
}

void ttkMultiTrackSource::addOneTrack(vtkUnstructuredGrid *dataset,
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
                                      vtkIntArray *pointTypeScalars) {
  std::vector<double> coords[3];

  for(size_t iDim = 0; iDim < 3; ++iDim) {
    auto &tab = coords[iDim];
    tab.resize(stopTimeExcl - startTime + NbSmooth);
    std::uniform_real_distribution<double> distrib(
      center[iDim] - 2 * baseStdDev[iDim], center[iDim] + 2 * baseStdDev[iDim]);
    for(auto &x : tab)
      x = distrib(generator);
    for(size_t iSmooth = 0; iSmooth < NbSmooth; ++iSmooth) {
      double last = tab.back();
      for(size_t jTab = tab.size(); jTab > 0; --jTab) {
        double newVal = (last + tab[jTab - 1]) / 2.;
        last = tab[jTab - 1];
        tab[jTab - 1] = newVal;
      }
      tab.pop_back();
    }
  }
  for(size_t jZ = 0; jZ < coords[2].size(); ++jZ)
    coords[2][jZ] += (startTime + jZ) * ZTranslation;

  bool isFirstPoint = true;
  auto insertPoint = [&](size_t iPt, double persVal) {
    points->InsertNextPoint(coords[0][iPt], coords[1][iPt], coords[2][iPt]);
    pointTypeScalars->InsertNextTuple1(
      (double)(int)ttk::CriticalType::Local_maximum);
    timeScalars->InsertNextTuple1(iPt + startTime);
    persistenceScalars->InsertNextTuple1(persVal);
    componentIds->InsertNextTuple1(componentId);
    valueScalars->InsertNextTuple1(0.);
    if(!isFirstPoint) {
      // Abusing permissivity of C++ is fun
      vtkIdType pts[2] = {points->GetNumberOfPoints() - 2, pts[0] + 1};
      dataset->InsertNextCell(VTK_LINE, 2, pts);
      costScalars->InsertNextTuple1(-1.);
      matchingIdScalars->InsertNextTuple1(dataset->GetNumberOfCells() - 1);
      lengthScalars->InsertNextTuple1(-1.); // Don't be relou
    }
    isFirstPoint = false;
    return;
  };
  size_t nbPoint = stopTimeExcl - startTime;

  /* Profile will look like this :
   *    ..
   *   /  \
   *  /    \
   * /      \
   *
   * With the central value repetated once or twice.
   * Central value is peakPers.
   */
  double increment
    = peakPers / ((nbPoint - 1) / 2); // Second division is Euclidean !
  for(size_t iPt = 0; iPt < nbPoint / 2; ++iPt)
    insertPoint(iPt, increment * iPt);
  for(size_t iPt = nbPoint / 2; iPt < nbPoint; ++iPt)
    insertPoint(iPt, increment * (nbPoint - 1 - iPt));

  return;
}

template <class URNG>
std::vector<size_t> randomMultiSelectFromRange(size_t k, size_t N, URNG &&g) {
  std::vector<size_t> ret;
  for(size_t i = 0; i < k; ++i)
    ret.push_back(std::uniform_int_distribution<size_t>(0, N)(g));

  std::sort(ret.begin(), ret.end());
  return ret;
}

int ttkMultiTrackSource::RequestData(vtkInformation *request,
                                     vtkInformationVector **inputVector,
                                     vtkInformationVector *outputVector) {
  if(outputMesh_)
    outputMesh_->Delete();
  outputMesh_ = vtkUnstructuredGrid::New();

  auto points = vtkSmartPointer<vtkPoints>::New();
  auto output = vtkUnstructuredGrid::GetData(outputVector);
  auto myOutput = vtkSmartPointer<vtkUnstructuredGrid>::New();
  myOutput->GetCellData()->Allocate(1);

  vtkSmartPointer<vtkDoubleArray> costScalars
    = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> persistenceScalars
    = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkDoubleArray> valueScalars
    = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkIntArray> matchingIdScalars
    = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> lengthScalars
    = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> timeScalars
    = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> componentIds
    = vtkSmartPointer<vtkIntArray>::New();
  vtkSmartPointer<vtkIntArray> pointTypeScalars
    = vtkSmartPointer<vtkIntArray>::New();
  costScalars->SetName("Cost");
  persistenceScalars->SetName("Persistence");
  valueScalars->SetName("Scalar");
  matchingIdScalars->SetName("MatchingIdentifier");
  lengthScalars->SetName("ComponentLength");
  timeScalars->SetName("TimeStep");
  componentIds->SetName("ConnectedComponentId");
  pointTypeScalars->SetName("CriticalType");

  double stddev[3] = {StdDev, StdDev, 0.};

  /** Old, parallel version
  std::uniform_int_distribution<> startDistrib(
    0, TotalLength - MinLength - (int)LengthPoisson);
  std::poisson_distribution<> poisson(LengthPoisson);
  for (size_t iLoc = 0; iLoc < NumberOfLocations; ++iLoc) {
    double phase = 2. * M_PI * iLoc / (double)NumberOfLocations;
    double center[3] = {Radius * std::cos(phase), Radius * std::sin(phase), 0.};
    for(size_t kTrack = 0; kTrack < NbTrackPerLocation; ++kTrack) {
      int start = startDistrib(generator);
      int end = start + MinLength + poisson(generator);
      addOneTrack(myOutput, points, center, stddev, start, end, 1.,
                  iLoc * NbTrackPerLocation + kTrack, costScalars,
                  persistenceScalars, valueScalars, matchingIdScalars,
                  lengthScalars, timeScalars, componentIds, pointTypeScalars);
    }
  }
  // */
  size_t massOfTracks = ceil(MassPortion * TotalLength);
  size_t massOfVoids = TotalLength - massOfTracks;
  if(massOfTracks < NbTrackPerLocation * MinLength) {
    printErr(
      "Parameter error : not enough mass of track to fulfill constraints");
    return 0;
  }
  massOfTracks -= NbTrackPerLocation * MinLength;
  for(size_t iLoc = 0; iLoc < NumberOfLocations; ++iLoc) {
    auto vecTrackSplits = randomMultiSelectFromRange(
      NbTrackPerLocation - 1, massOfTracks, generator);
    vecTrackSplits.insert(vecTrackSplits.begin(), 0);
    vecTrackSplits.push_back(massOfTracks);

    auto vecVoidSplits
      = randomMultiSelectFromRange(NbTrackPerLocation, massOfVoids, generator);

    double phase = 2. * M_PI * iLoc / (double)NumberOfLocations;
    double center[3] = {Radius * std::cos(phase), Radius * std::sin(phase), 0.};

    for(size_t kTrack = 0; kTrack < NbTrackPerLocation; ++kTrack) {
      int start
        = vecVoidSplits[kTrack] + vecTrackSplits[kTrack] + kTrack * MinLength;
      int end = vecVoidSplits[kTrack] + vecTrackSplits[kTrack + 1]
                + (kTrack + 1) * MinLength;
      addOneTrack(myOutput, points, center, stddev, start, end, 1.,
                  iLoc * NbTrackPerLocation + kTrack, costScalars,
                  persistenceScalars, valueScalars, matchingIdScalars,
                  lengthScalars, timeScalars, componentIds, pointTypeScalars);
    }
  }

  myOutput->SetPoints(points);
  myOutput->GetCellData()->AddArray(costScalars);
  myOutput->GetPointData()->AddArray(persistenceScalars);
  myOutput->GetPointData()->AddArray(valueScalars);
  myOutput->GetCellData()->AddArray(matchingIdScalars);
  myOutput->GetCellData()->AddArray(lengthScalars);
  myOutput->GetPointData()->AddArray(timeScalars);
  myOutput->GetPointData()->AddArray(componentIds);
  myOutput->GetPointData()->AddArray(pointTypeScalars);

  outputMesh_->ShallowCopy(myOutput);
  output->ShallowCopy(outputMesh_);
  return 1;
}
