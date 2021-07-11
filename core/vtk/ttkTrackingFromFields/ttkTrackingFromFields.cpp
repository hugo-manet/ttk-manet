#include "VineyardTracking.hpp"
#include <vtkInformation.h>

#include <ttkMacros.h>
#include <ttkTrackingFromFields.h>
#include <ttkTrackingFromPersistenceDiagrams.h>
#include <ttkUtils.h>

vtkStandardNewMacro(ttkTrackingFromFields);

constexpr unsigned long long str2int(const char *str, int h = 0) {
  return !str[h] ? 5381 : (str2int(str, h + 1) * 33) ^ str[h];
}

ttkTrackingFromFields::ttkTrackingFromFields() {
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkTrackingFromFields::FillOutputPortInformation(int port,
                                                     vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}
int ttkTrackingFromFields::FillInputPortInformation(int port,
                                                    vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
    return 1;
  }
  return 0;
}

// (*) Vineyard-driven approach
template <class dataType, class triangulationType>
int ttkTrackingFromFields::trackWithVineyards(
  vtkDataSet *input,
  vtkUnstructuredGrid *output,
  unsigned long fieldNumber,
  triangulationType *triangulation) {

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkUnstructuredGrid> persistenceDiagram
    = vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkSmartPointer<vtkDoubleArray> costScalars
    = vtkSmartPointer<vtkDoubleArray>::New();
  vtkSmartPointer<vtkIntArray> matchTypeScalars
    = vtkSmartPointer<vtkIntArray>::New();
  costScalars->SetName("Cost");
  matchTypeScalars->SetName("MatchType");

  SimplexId nbNodes = triangulation->getNumberOfVertices();
  triangulation->preconditionVertexNeighbors();

#ifndef NONOISE
  srand(42); // Try 42 for bug at time 0.00996
  for(int idFieldStart = 0; idFieldStart < fieldNumber; ++idFieldStart)
    for(int i = 0; i < nbNodes; ++i)
      ((double *)inputData_[idFieldStart])[i]
        += (rand() - RAND_MAX / 2.) / (RAND_MAX * 1000.);
#endif

  double globalCost = 0.0;
#pragma omp parallel for
  for(int idFieldStart = 0; idFieldStart < fieldNumber - 1; ++idFieldStart) {
    auto res = ttk::buildTree(triangulation, (double *)inputData_[idFieldStart],
                              (double *)inputData_[idFieldStart + 1]);
    auto &nodesVec = res.first;
    ttk::EventQueue &events = res.second;

    std::map<SimplexId, std::pair<double, double>> startPairs;
    for(SimplexId i = 0; i < nodesVec.size(); ++i) {
      if(nodesVec[i].pairOfMax != NULL)
        startPairs.insert({i,
                           {nodesVec[i].pairOfMax->saddle->scalarStart,
                            nodesVec[i].pairOfMax->max->scalarEnd}});
    }

    ttk::loopQueue(events);
    double localCost = 0.0;

#pragma omp critical
    {
      auto addPair = [&](SimplexId start, SimplexId end,
                         std::pair<double, double> startPair,
                         std::pair<double, double> endPair) {
        int typeOfMatch = 1; // NORMAL
        if(start == -1) {
          start = end;
          typeOfMatch = 0; // CREATION
        }
        if(end == -1) {
          end = start;
          typeOfMatch = 2; // DELETION
        }
        float startCoord[3], endCoord[3];
        triangulation->getVertexPoint(
          start, startCoord[0], startCoord[1], startCoord[2]);
        triangulation->getVertexPoint(
          end, endCoord[0], endCoord[1], endCoord[2]);

        if(UseGeometricSpacing) {
          startCoord[2] += Spacing * (idFieldStart);
          endCoord[2] += Spacing * (idFieldStart + 1);
        }
        points->InsertNextPoint(startCoord[0], startCoord[1], startCoord[2]);
        points->InsertNextPoint(endCoord[0], endCoord[1], endCoord[2]);
        vtkIdType line[2]
          = {points->GetNumberOfPoints() - 2, points->GetNumberOfPoints() - 1};
        persistenceDiagram->InsertNextCell(VTK_LINE, 2, line);
        matchTypeScalars->InsertNextTuple1(typeOfMatch);
        double costSquared = (startPair.first - endPair.first)
                               * (startPair.first - endPair.first)
                             + (startPair.second - endPair.second)
                                 * (startPair.second - endPair.second);
        localCost += costSquared;
        globalCost += costSquared;
        costScalars->InsertNextTuple1(std::sqrt(costSquared));
      };
      for(SimplexId i = 0; i < nodesVec.size(); ++i) {
        if(nodesVec[i].pairOfMax != NULL) {
          std::pair<double, double> endPair{
            nodesVec[i].pairOfMax->saddle->scalarEnd,
            nodesVec[i].pairOfMax->max->scalarEnd};
          double middle = (endPair.first + endPair.second) / 2.;
          std::pair<double, double> startPair{middle, middle};
          if(nodesVec[i].pairOfMax->idFirstMax != -1) {
            auto it = startPairs.find(nodesVec[i].pairOfMax->idFirstMax);
            startPair = it->second;
            startPairs.erase(it);
          }
          addPair(nodesVec[i].pairOfMax->idFirstMax, i, startPair, endPair);
        }
      }
      for(auto it : startPairs) {
        auto stp = it.second;
        double middle = (stp.first + stp.second) / 2.;
        std::pair<double, double> endPair{middle, middle};
        addPair(-1, it.first, stp, endPair);
      }
      std::cout << "Local cost for match " << idFieldStart << " : "
                << globalCost << std::endl;
    }
  }

  std::cout << "Global cost : " << globalCost << std::endl;

  persistenceDiagram->SetPoints(points);
  persistenceDiagram->GetCellData()->AddArray(costScalars);
  persistenceDiagram->GetCellData()->AddArray(matchTypeScalars);

  output->ShallowCopy(persistenceDiagram);
  return 1;
}

// (*) Persistence-driven approach
template <class dataType, class triangulationType>
int ttkTrackingFromFields::trackWithPersistenceMatching(
  vtkDataSet *input,
  vtkUnstructuredGrid *output,
  unsigned long fieldNumber,
  const triangulationType *triangulation) {

  using trackingTuple = ttk::trackingTuple;

  // 1. get persistence diagrams.
  std::vector<std::vector<diagramTuple>> persistenceDiagrams(
    fieldNumber, std::vector<diagramTuple>());

  this->performDiagramComputation<dataType, triangulationType>(
    (int)fieldNumber, persistenceDiagrams, triangulation);

  // 2. call feature tracking with threshold.
  std::vector<std::vector<matchingTuple>> outputMatchings(
    fieldNumber - 1, std::vector<matchingTuple>());

  double spacing = Spacing;
  std::string algorithm = DistanceAlgorithm;
  double alpha = Alpha;
  double tolerance = Tolerance;
  bool is3D = true; // Is3D;
  std::string wasserstein = WassersteinMetric;

  ttk::TrackingFromPersistenceDiagrams tfp{};
  tfp.setThreadNumber(this->threadNumber_);
  tfp.performMatchings(
    (int)fieldNumber, persistenceDiagrams, outputMatchings,
    algorithm, // Not from paraview, from enclosing tracking plugin
    wasserstein, tolerance, is3D,
    alpha, // Blending
    PX, PY, PZ, PS, PE // Coefficients
  );

  vtkNew<vtkPoints> points{};
  vtkNew<vtkUnstructuredGrid> persistenceDiagram{};

  vtkNew<vtkDoubleArray> costScalars{};
  vtkNew<vtkDoubleArray> persistenceScalars{};
  vtkNew<vtkDoubleArray> valueScalars{};
  vtkNew<vtkIntArray> matchingIdScalars{};
  vtkNew<vtkIntArray> lengthScalars{};
  vtkNew<vtkIntArray> timeScalars{};
  vtkNew<vtkIntArray> componentIds{};
  vtkNew<vtkIntArray> pointTypeScalars{};

  costScalars->SetName("Cost");
  persistenceScalars->SetName("Persistence");
  valueScalars->SetName("Scalar");
  matchingIdScalars->SetName("MatchingIdentifier");
  lengthScalars->SetName("ComponentLength");
  timeScalars->SetName("TimeStep");
  componentIds->SetName("ConnectedComponentId");
  pointTypeScalars->SetName("CriticalType");

  // (+ vertex id)
  std::vector<trackingTuple> trackingsBase;
  tfp.performTracking(persistenceDiagrams, outputMatchings, trackingsBase);

  std::vector<std::set<int>> trackingTupleToMerged(
    trackingsBase.size(), std::set<int>());

  if(DoPostProc) {
    tfp.performPostProcess(persistenceDiagrams, trackingsBase,
                           trackingTupleToMerged, PostProcThresh);
  }

  bool useGeometricSpacing = UseGeometricSpacing;

  // Build mesh.
  ttkTrackingFromPersistenceDiagrams::buildMesh(
    trackingsBase, outputMatchings, persistenceDiagrams, useGeometricSpacing,
    spacing, DoPostProc, trackingTupleToMerged, points, persistenceDiagram,
    costScalars, persistenceScalars, valueScalars, matchingIdScalars,
    lengthScalars, timeScalars, componentIds, pointTypeScalars);

  output->ShallowCopy(persistenceDiagram);

  return 1;
}

int ttkTrackingFromFields::RequestData(vtkInformation *request,
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {

  auto input = vtkDataSet::GetData(inputVector[0]);
  auto output = vtkUnstructuredGrid::GetData(outputVector);

  ttk::Triangulation *triangulation = ttkAlgorithm::GetTriangulation(input);
  if(!triangulation)
    return 0;

  this->preconditionTriangulation(triangulation);

  // Test validity of datasets
  if(input == nullptr || output == nullptr) {
    return -1;
  }

  // Get number and list of inputs.
  std::vector<vtkDataArray *> inputScalarFieldsRaw;
  std::vector<vtkDataArray *> inputScalarFields;
  const auto pointData = input->GetPointData();
  int numberOfInputFields = pointData->GetNumberOfArrays();
  if(numberOfInputFields < 3) {
    this->printErr("Not enough input fields to perform tracking.");
  }

  vtkDataArray *firstScalarField = pointData->GetArray(0);

  for(int i = 0; i < numberOfInputFields; ++i) {
    vtkDataArray *currentScalarField = pointData->GetArray(i);
    if(currentScalarField == nullptr
       || currentScalarField->GetName() == nullptr) {
      continue;
    }
    std::string sfname{currentScalarField->GetName()};
    if(sfname.rfind("_Order") == (sfname.size() - 6)) {
      continue;
    }
    if(firstScalarField->GetDataType() != currentScalarField->GetDataType()) {
      this->printErr("Inconsistent field data type or size between fields `"
                     + std::string{firstScalarField->GetName()} + "' and `"
                     + sfname + "'");
      return -1;
    }
    inputScalarFieldsRaw.push_back(currentScalarField);
  }

  std::sort(inputScalarFieldsRaw.begin(), inputScalarFieldsRaw.end(),
            [](vtkDataArray *a, vtkDataArray *b) {
              std::string s1 = a->GetName();
              std::string s2 = b->GetName();
              return std::lexicographical_compare(
                s1.begin(), s1.end(), s2.begin(), s2.end());
            });

  numberOfInputFields = inputScalarFieldsRaw.size();
  int end = EndTimestep <= 0 ? numberOfInputFields
                             : std::min(numberOfInputFields, EndTimestep);
  for(int i = StartTimestep; i < end; i += Sampling) {
    vtkDataArray *currentScalarField = inputScalarFieldsRaw[i];
    // Print scalar field names:
    // std::cout << currentScalarField->GetName() << std::endl;
    inputScalarFields.push_back(currentScalarField);
  }

  // Input -> persistence filter.
  std::string algorithm = DistanceAlgorithm;
  int pvalg = PVAlgorithm;
  bool useTTKMethod = false;

  if(pvalg >= 0) {
    switch(pvalg) {
      case 0:
      case 1:
      case 2:
      case 3:
        useTTKMethod = true;
        break;
      case 4:
      case 5:
        break;
      default:
        this->printMsg("Unrecognized tracking method.");
        break;
    }
  } else {
    switch(str2int(algorithm.c_str())) {
      case str2int("0"):
      case str2int("ttk"):
      case str2int("1"):
      case str2int("legacy"):
      case str2int("2"):
      case str2int("geometric"):
      case str2int("3"):
      case str2int("parallel"):
        useTTKMethod = true;
        break;
      case str2int("4"):
      case str2int("greedy"):
      case str2int("5"):
      case str2int("vineyard"):
        break;
      default:
        this->printMsg("Unrecognized tracking method.");
        break;
    }
  }

  // 0. get data
  int fieldNumber = inputScalarFields.size();
  std::vector<void *> inputFields(fieldNumber);
  for(int i = 0; i < fieldNumber; ++i) {
    inputFields[i] = ttkUtils::GetVoidPointer(inputScalarFields[i]);
  }
  this->setInputScalars(inputFields);

  // 0'. get offsets
  std::vector<ttk::SimplexId *> inputOrders(fieldNumber);
  for(int i = 0; i < fieldNumber; ++i) {
    this->SetInputArrayToProcess(0, 0, 0, 0, inputScalarFields[i]->GetName());
    auto orderArray = this->GetOrderArray(input, 0, 0, false);
    inputOrders[i]
      = static_cast<ttk::SimplexId *>(ttkUtils::GetVoidPointer(orderArray));
  }
  this->setInputOffsets(inputOrders);

  int status = 0;
  if(useTTKMethod) {
    ttkVtkTemplateMacro(
      inputScalarFields[0]->GetDataType(), triangulation->getType(),
      (status = this->trackWithPersistenceMatching<VTK_TT, TTK_TT>(
         input, output, fieldNumber, (TTK_TT *)triangulation->getData())));
  } else {
    this->printMsg("Experimental tracking method...");
    ttkVtkTemplateMacro(
      inputScalarFields[0]->GetDataType(), triangulation->getType(),
      (status = this->trackWithVineyards<VTK_TT, TTK_TT>(
         input, output, fieldNumber, (TTK_TT *)triangulation->getData())));
  }

  return status;
}
