#include <map>

#include "DataTypes.h"
#include "TracksMatching.h"
#include <sstream>
#include <ttkTracksMatching.h>

#include <vtkInformation.h>

#include <vtkCellData.h>
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
    double realZ = points->GetPoint(i)[2];
    if(ZTranslation < 0 && timeStep != 0.)
      ZTranslation = realZ / timeStep;
    realZ -= ZTranslation * timeStep;

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
  auto outputDataSet = vtkUnstructuredGrid::GetData(outputVector, 0);
  outputDataSet->AllocateEstimate(
    inputTracksBidder->GetNumberOfCells() + inputTracksGood->GetNumberOfCells(),
    2);
  auto outputPoints = vtkSmartPointer<vtkPoints>::New();
  outputDataSet->SetPoints(outputPoints);
  vtkIdType nbPoints = 0;

  this->bidders = getTracksFromObject(inputTracksBidder);
  this->goods = getTracksFromObject(inputTracksGood);

  const double p = this->PowerParameter;
  const double lambda_p = std::pow(this->LambdaParameter, p);
  const double phi_p = std::pow(0.5, p - 1.);
  const double tN = TimeNormalization;
  const double gL = GeometricalLifting;
  this->run(this->PowerParameter, lambda_p, tN, gL);

  auto deletionTypeOutputArray = vtkSmartPointer<vtkIntArray>::New();
  enum deletionType { MATCHED, COMPRESSED, DELETED };
  deletionTypeOutputArray->SetName("DeletionType");
  outputDataSet->GetCellData()->AddArray(deletionTypeOutputArray);
  auto localCostOutputArray = vtkSmartPointer<vtkDoubleArray>::New();
  localCostOutputArray->SetName("LocalCost");
  outputDataSet->GetCellData()->AddArray(localCostOutputArray);
  auto trackCostOutputArray = vtkSmartPointer<vtkDoubleArray>::New();
  trackCostOutputArray->SetName("TrackCost");
  outputDataSet->GetCellData()->AddArray(trackCostOutputArray);

  auto insertPoint = [&](TimedPoint pt, ttk::CriticalType ct, bool offset) {
    outputPoints->InsertNextPoint(pt.x, pt.y + offset * this->DisplayOffset,
                                  pt.z + pt.timeStep * this->ZTranslation);

    ++nbPoints;
  };

  auto addTrackDeleteToOutput = [&](const Track &track, bool offset) {
    if(track.empty())
      return 0.;
    insertPoint(track[0], track.trackType, offset);
    double trackDelCost = track.deletionCost(p);
    for(size_t iPt = 1; iPt < track.size(); ++iPt) {
      insertPoint(track[iPt], track.trackType, offset);
      vtkIdType deletedTrackSegment[2] = {nbPoints - 2, nbPoints - 1};

      outputDataSet->InsertNextCell(VTK_LINE, 2, deletedTrackSegment);
      deletionTypeOutputArray->InsertNextValue(DELETED);
      localCostOutputArray->InsertNextValue(std::pow(
        ttk::sumPow(p, track[iPt - 1].persistence, track[iPt].persistence),
        1. / p));
      trackCostOutputArray->InsertNextValue(trackDelCost);
    }
    return trackDelCost;
  };

  auto addTrackMatchToOutput = [&](TWED &match) {
    ttk::CriticalType trackType = match.bidderTrack->trackType;
    vtkIdType iPtBid = nbPoints;
    for(auto &pt : *match.bidderTrack)
      insertPoint(pt, trackType, false);
    vtkIdType jPtGood = nbPoints;
    for(auto &pt : *match.goodTrack)
      insertPoint(pt, trackType, true);

    auto localDelCost = [&](Track &track, size_t i) {
      return std::pow(
        lambda_p + phi_p * ttk::distPow(track[i], track[i + 1], p, tN, gL),
        1. / p);
    };
    for(auto [iB, jG] : match.matchedSegments)
      if(iB == -1) {
        vtkIdType compressedTrackSegment[2] = {jPtGood + jG, jPtGood + jG + 1};
        outputDataSet->InsertNextCell(VTK_LINE, 2, compressedTrackSegment);
        deletionTypeOutputArray->InsertNextValue(COMPRESSED);
        localCostOutputArray->InsertNextValue(
          localDelCost(*match.goodTrack, jG));
        trackCostOutputArray->InsertNextValue(match.value);
      } else if(jG == -1) {
        vtkIdType compressedTrackSegment[2] = {iPtBid + iB, iPtBid + iB + 1};
        outputDataSet->InsertNextCell(VTK_LINE, 2, compressedTrackSegment);
        deletionTypeOutputArray->InsertNextValue(COMPRESSED);
        localCostOutputArray->InsertNextValue(
          localDelCost(*match.bidderTrack, iB));
        trackCostOutputArray->InsertNextValue(match.value);
      } else {
        vtkIdType matchedTrackSegment[4]
          = {iPtBid + iB, iPtBid + iB + 1, jPtGood + jG + 1, jPtGood + jG};
        outputDataSet->InsertNextCell(VTK_QUAD, 4, matchedTrackSegment);
        deletionTypeOutputArray->InsertNextValue(MATCHED);
        localCostOutputArray->InsertNextValue(
          std::pow(ttk::distPow((*match.bidderTrack)[iB],
                                (*match.goodTrack)[jG], p, tN, gL)
                     + ttk::distPow((*match.bidderTrack)[iB + 1],
                                    (*match.goodTrack)[jG + 1], p, tN, gL),
                   1. / p));
        trackCostOutputArray->InsertNextValue(match.value);
      }
    return match.value;
  };

  double globalCost = 0.;
  for(auto [iB, jG] : matchedTracks) {
    if(iB == -1)
      globalCost += std::pow(addTrackDeleteToOutput(goods[jG], true),p);
    else if(jG == -1)
      globalCost += std::pow(addTrackDeleteToOutput(bidders[iB], false),p);
    else
      globalCost += std::pow(addTrackMatchToOutput(distances(iB, jG)),p);
  }

  std::stringstream globCostStringStream;
  globCostStringStream << "Global cost : " << std::pow(globalCost, 1. / p);

  printMsg(globCostStringStream.str());

  // return success
  return 1;
}
