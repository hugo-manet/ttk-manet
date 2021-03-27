#include <map>

#include "DataTypes.h"
#include "TracksMatching.h"
#include <ttkTracksMatching.h>

#include <vtkInformation.h>

#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>

#include <ttkMacros.h>
#include <ttkUtils.h>
#include <vtkUnstructuredGrid.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkTracksMatching);

/**
 * TODO 7: Implement the filter constructor and destructor in the cpp file.
 *
 * The constructor has to specify the number of input and output ports
 * with the functions SetNumberOfInputPorts and SetNumberOfOutputPorts,
 * respectively. It should also set default values for all filter
 * parameters.
 *
 * The destructor is usually empty unless you want to manage memory
 * explicitly, by for example allocating memory on the heap that needs
 * to be freed when the filter is destroyed.
 */
ttkTracksMatching::ttkTracksMatching() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

ttkTracksMatching::~ttkTracksMatching() {
}

/**
 * TODO 8: Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkTracksMatching::FillInputPortInformation(int port,
                                                vtkInformation *info) {
  if(port == 0 || port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

/**
 * TODO 9: Specify the data object type of each output port
 *
 * This method specifies in the port information object the data type of the
 * corresponding output objects. It is possible to either explicitly
 * specify a type by adding a vtkDataObject::DATA_TYPE_NAME() key:
 *
 *      info->Set( vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid" );
 *
 * or to pass a type of an input port to an output port by adding the
 * ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT() key (see below).
 *
 * Note: prior to the execution of the RequestData method the pipeline will
 * initialize empty output data objects based on this information.
 */
int ttkTracksMatching::FillOutputPortInformation(int port,
                                                 vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

std::vector<ttk::TracksMatching::Track>
  ttkTracksMatching::getTracksFromObject(vtkUnstructuredGrid *obj) {
  using Track = ttk::TracksMatching::Track;
  using TimedPoint = ttk::TracksMatching::TimedPoint;
  std::vector<Track> ret;
  /** doesn't exist ?
  ttkSimplexIdTypeArray *vertexIdentifierScalars
    = ttkSimplexIdTypeArray::SafeDownCast(
      obj->GetPointData()->GetArray(ttk::VertexScalarFieldName));
  // */
  vtkIntArray *trackIdScalars = vtkIntArray::SafeDownCast(
    obj->GetPointData()->GetArray("ConnectedComponentId"));

  vtkDoubleArray *persistenceScalars = vtkDoubleArray::SafeDownCast(
    obj->GetPointData()->GetArray("Persistence"));

  vtkIntArray *timeStepScalars
    = vtkIntArray::SafeDownCast(obj->GetPointData()->GetArray("TimeStep"));

  vtkDoubleArray *valueScalars
    = vtkDoubleArray::SafeDownCast(obj->GetPointData()->GetArray("Scalar"));

  vtkIntArray *criticalTypeScalars
    = vtkIntArray::SafeDownCast(obj->GetPointData()->GetArray("CriticalType"));

  if(!trackIdScalars || !persistenceScalars || !timeStepScalars || !valueScalars
     || !criticalTypeScalars)
    throw "No something";

  vtkPoints *points = obj->GetPoints();
  size_t nbOfPoints = obj->GetNumberOfPoints();

  if(nbOfPoints == 0)
    throw 0;

  auto pointFromIndex = [&](size_t i) -> TimedPoint {
    double timeStep = timeStepScalars->GetValue(i);
    double realZ = (ZTranslation < 0)
                     ? 0
                     : points->GetPoint(i)[2] - ZTranslation * timeStep;
    return TimedPoint(timeStep, valueScalars->GetValue(i),
                      points->GetPoint(i)[0], points->GetPoint(i)[1], realZ,
                      persistenceScalars->GetValue(i));
  };
  std::map<int, int> trackId; // Tracks could be filtered by user
  int nextTrackId = 0;
  for(size_t i = 0; i < nbOfPoints; ++i) {
    int tID = trackIdScalars->GetValue(i);
    if(trackId.find(tID) == trackId.end()) {
      trackId[tID] = nextTrackId++;
      ret.push_back(Track());
      ret.back().trackType
        = (ttk::CriticalType)criticalTypeScalars->GetValue(i);
    }
    ret[trackId[tID]].push_back(pointFromIndex(i));
  }
  for(auto &track : ret) {
    std::sort(track.begin(), track.end(), [](TimedPoint &tA, TimedPoint &tB) {
      return tA.timeStep < tB.timeStep;
    });
    // remove duplicate from tracks
    size_t iWriten = 0;
    for(size_t iRead = 1; iRead < track.size(); ++iRead)
      if(track[iRead].timeStep != track[iWriten].timeStep)
        track[++iWriten] = track[iRead];
    track.erase(track.begin() + iWriten + 1, track.end());
  }

  return ret;
}

int ttkTracksMatching::RequestData(vtkInformation *request,
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector) {

  // Get input object from input vector
  // Note: has to be a vtkUnstrcturedGrid as required by
  // FillInputPortInformation
  vtkUnstructuredGrid *inputTracksBidder
    = vtkUnstructuredGrid::GetData(inputVector[0]);
  vtkUnstructuredGrid *inputTracksGood
    = vtkUnstructuredGrid::GetData(inputVector[1]);
  if(!inputTracksBidder || !inputTracksGood)
    return 0;

  // Get output vtkDataSet (which was already instantiated based on the
  // information provided by FillOutputPortInformation)
  vtkDataSet *outputDataSet = vtkUnstructuredGrid::GetData(outputVector, 0);

  this->bidders = getTracksFromObject(inputTracksBidder);
  this->goods = getTracksFromObject(inputTracksGood);

  this->run(this->PowerParameter, 1, 1, 1);

  // return success
  return 1;
}
