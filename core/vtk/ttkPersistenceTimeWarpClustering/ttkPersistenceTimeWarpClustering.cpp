#include <ttkMacros.h>
#include <ttkPersistenceTimeWarpClustering.h>
#include <ttkUtils.h>
#include <vtkFieldData.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPersistenceTimeWarpClustering)

  ttkPersistenceTimeWarpClustering::ttkPersistenceTimeWarpClustering() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(3);
}

int ttkPersistenceTimeWarpClustering::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
    info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), 1);
  } else
    return 0;
  return 1;
}

int ttkPersistenceTimeWarpClustering::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0 || port == 1)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
  else if(port == 2)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  else
    return 0;
  return 1;
}

// to adapt if your wrapper does not inherit from vtkDataSetAlgorithm
int ttkPersistenceTimeWarpClustering::RequestData(
  vtkInformation * /*request*/,
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  Memory m;

  // Get input data

  auto nCurves = inputVector[0]->GetNumberOfInformationObjects();
  std::vector<vtkMultiBlockDataSet *> blocks(nCurves);
  std::vector<std::vector<vtkUnstructuredGrid *>> inputDiagramGrids(nCurves);

  if(nCurves == 0) {
    this->printErr("Input is empty");
    return 0;
  }

  // number of diagrams per input block
  std::vector<size_t> nDiagOfCurve(nCurves);

  for(int iCurve = 0; iCurve < nCurves; ++iCurve) {
    blocks[iCurve] = vtkMultiBlockDataSet::GetData(inputVector[0], iCurve);
    if(blocks[iCurve] != nullptr) {
      nDiagOfCurve[iCurve] = blocks[iCurve]->GetNumberOfBlocks();
      for(size_t jDiag = 0; jDiag < nDiagOfCurve[iCurve]; ++jDiag) {
        inputDiagramGrids[iCurve].emplace_back(
          vtkUnstructuredGrid::SafeDownCast(blocks[iCurve]->GetBlock(jDiag)));
      }
    }
  }

  // Sanity check
  for(const auto &curveGrid : inputDiagramGrids) {
    if(nDiagOfCurve[0] != curveGrid.size()) {
      this->printErr("Input curves aren't all the same size. Fatal for now");
      return 0;
    }
    for(const auto &vtu : curveGrid)
      if(vtu == nullptr) {
        this->printErr("Input diagrams are not all vtkUnstructuredGrid");
        return 0;
      }
  }

  std::vector<DiagramCurve> inputCurves;
  for(const auto &curveGrid : inputDiagramGrids) {
    inputCurves.emplace_back();
    for(const auto &vtu : curveGrid) {
      inputCurves.back().emplace_back();
      this->getPersistenceDiagram(inputCurves.back().back(), vtu);
    }
  }
  DiagramCurve barycenter(nDiagOfCurve[0]);
  using dataType = double;
  std::vector<std::vector<std::vector<matchingTuple>>> all_matchings_;
  this->execute(inputCurves, barycenter, all_matchings_);

  // Set outputs
  auto outputInitialDiagrams = vtkMultiBlockDataSet::SafeDownCast(
    outputVector->GetInformationObject(0)->Get(vtkDataObject::DATA_OBJECT()));
  auto outputBarycenterCurves = vtkMultiBlockDataSet::SafeDownCast(
    outputVector->GetInformationObject(1)->Get(vtkDataObject::DATA_OBJECT()));
  auto outputMatching = vtkUnstructuredGrid::SafeDownCast(
    outputVector->GetInformationObject(2)->Get(vtkDataObject::DATA_OBJECT()));

  // Copy input to output, add curve and diagram index, and store diagram
  std::vector<int> firstDiagramIDOfCurve;
  size_t indexOfDiag = 0;
  for(int iCurve = 0; iCurve < nCurves; ++iCurve) {
    firstDiagramIDOfCurve.push_back(indexOfDiag);
    auto curveBlock = vtkSmartPointer<vtkMultiBlockDataSet>::New();
    for(int jDiag = 0; jDiag < nDiagOfCurve[iCurve]; ++jDiag) {
      auto &diag = inputDiagramGrids[iCurve][jDiag];
      auto copiedDiag
        = vtkSmartPointer<vtkUnstructuredGrid>::Take(diag->NewInstance());
      copiedDiag->ShallowCopy(diag);

      vtkNew<vtkIntArray> curveIndex_p{};
      curveIndex_p->SetName("CurveID");
      curveIndex_p->SetNumberOfTuples(copiedDiag->GetNumberOfPoints());
      curveIndex_p->FillValue(iCurve);
      copiedDiag->GetPointData()->AddArray(curveIndex_p);
      vtkNew<vtkIntArray> curveIndex_c{};
      curveIndex_c->SetName("CurveID");
      curveIndex_c->SetNumberOfTuples(copiedDiag->GetNumberOfCells());
      curveIndex_c->FillValue(iCurve);
      copiedDiag->GetCellData()->AddArray(curveIndex_c);

      vtkNew<vtkIntArray> diagramIndex_p{};
      diagramIndex_p->SetName("DiagramID");
      diagramIndex_p->SetNumberOfTuples(copiedDiag->GetNumberOfPoints());
      diagramIndex_p->FillValue(indexOfDiag);
      copiedDiag->GetPointData()->AddArray(diagramIndex_p);
      vtkNew<vtkIntArray> diagramIndex_c{};
      diagramIndex_c->SetName("DiagramID");
      diagramIndex_c->SetNumberOfTuples(copiedDiag->GetNumberOfCells());
      diagramIndex_c->FillValue(indexOfDiag);
      copiedDiag->GetCellData()->AddArray(diagramIndex_c);

      curveBlock->SetBlock(jDiag, copiedDiag);
      ++indexOfDiag;
    }
    outputInitialDiagrams->SetBlock(iCurve, curveBlock);
  }

  return 1;
}

double ttkPersistenceTimeWarpClustering::getPersistenceDiagram(
  ttk::Diagram &diagram, vtkUnstructuredGrid *CTPersistenceDiagram_) {
  vtkIntArray *vertexIdentifierScalars
    = vtkIntArray::SafeDownCast(CTPersistenceDiagram_->GetPointData()->GetArray(
      ttk::VertexScalarFieldName));

  vtkIntArray *nodeTypeScalars = vtkIntArray::SafeDownCast(
    CTPersistenceDiagram_->GetPointData()->GetArray("CriticalType"));

  vtkIntArray *pairIdentifierScalars = vtkIntArray::SafeDownCast(
    CTPersistenceDiagram_->GetCellData()->GetArray("PairIdentifier"));

  vtkIntArray *extremumIndexScalars = vtkIntArray::SafeDownCast(
    CTPersistenceDiagram_->GetCellData()->GetArray("PairType"));

  const auto persistenceScalars = vtkDoubleArray::SafeDownCast(
    CTPersistenceDiagram_->GetCellData()->GetArray("Persistence"));
  const auto birthScalars = vtkDoubleArray::SafeDownCast(
    CTPersistenceDiagram_->GetPointData()->GetArray("Birth"));
  const auto deathScalars = vtkDoubleArray::SafeDownCast(
    CTPersistenceDiagram_->GetPointData()->GetArray("Death"));

  const auto critCoordinates = vtkFloatArray::SafeDownCast(
    CTPersistenceDiagram_->GetPointData()->GetArray("Coordinates"));
  const auto points = CTPersistenceDiagram_->GetPoints();

  const bool embed = birthScalars != nullptr && deathScalars != nullptr;

#ifndef TTK_ENABLE_KAMIKAZE
  if(!embed && critCoordinates == nullptr) {
    // missing data
    return -2;
  }
#endif // TTK_ENABLE_KAMIKAZE

  int pairingsSize = (int)pairIdentifierScalars->GetNumberOfTuples();
  // FIX : no more missed pairs
  for(int pair_index = 0; pair_index < pairingsSize; pair_index++) {
    const float index_of_pair = pair_index;
    if(*pairIdentifierScalars->GetTuple(pair_index) != -1)
      pairIdentifierScalars->SetTuple(pair_index, &index_of_pair);
  }

#ifndef TTK_ENABLE_KAMIKAZE
  if(pairingsSize < 1 || !vertexIdentifierScalars || !pairIdentifierScalars
     || !nodeTypeScalars || !persistenceScalars || !extremumIndexScalars
     || !points) {
    return -2;
  }
#endif // TTK_ENABLE_KAMIKAZE

  if(NumberOfClusters == 1) {
    diagram.resize(pairingsSize);
  } else {
    diagram.resize(pairingsSize + 1);
  }
  int nbNonCompact = 0;
  double max_dimension = 0;

  // skip diagonal cell (corresponding points already dealt with)
  for(int i = 0; i < pairingsSize - 1; ++i) {

    int vertexId1 = vertexIdentifierScalars->GetValue(2 * i);
    int vertexId2 = vertexIdentifierScalars->GetValue(2 * i + 1);
    int nodeType1 = nodeTypeScalars->GetValue(2 * i);
    int nodeType2 = nodeTypeScalars->GetValue(2 * i + 1);

    int pairIdentifier = pairIdentifierScalars->GetValue(i);
    int pairType = extremumIndexScalars->GetValue(i);
    double persistence = persistenceScalars->GetValue(i);

    std::array<double, 3> coordsBirth{}, coordsDeath{};

    const auto i0 = 2 * i;
    const auto i1 = 2 * i + 1;

    double birth, death;

    if(embed) {
      points->GetPoint(i0, coordsBirth.data());
      points->GetPoint(i1, coordsDeath.data());
      birth = birthScalars->GetValue(i0);
      death = deathScalars->GetValue(i1);
    } else {
      critCoordinates->GetTuple(i0, coordsBirth.data());
      critCoordinates->GetTuple(i1, coordsDeath.data());
      birth = points->GetPoint(i0)[0];
      death = points->GetPoint(i1)[1];
    }

    if(pairIdentifier != -1 && pairIdentifier < pairingsSize) {
      if(pairIdentifier == 0) {
        max_dimension = persistence;

        if(NumberOfClusters == 1) {
          diagram[0] = std::make_tuple(
            vertexId1, CriticalType::Local_minimum, vertexId2,
            CriticalType::Local_maximum, persistence, pairType, birth,
            coordsBirth[0], coordsBirth[1], coordsBirth[2], death,
            coordsDeath[0], coordsDeath[1], coordsDeath[2]);
        } else {
          diagram[0] = std::make_tuple(
            vertexId1, CriticalType::Local_minimum, vertexId2,
            CriticalType::Saddle1, persistence, pairType, birth, coordsBirth[0],
            coordsBirth[1], coordsBirth[2], death, coordsDeath[0],
            coordsDeath[1], coordsDeath[2]);
          diagram[pairingsSize] = std::make_tuple(
            vertexId1, CriticalType::Saddle1, vertexId2,
            CriticalType::Local_maximum, persistence, pairType, birth,
            coordsBirth[0], coordsBirth[1], coordsBirth[2], death,
            coordsDeath[0], coordsDeath[1], coordsDeath[2]);
        }

      } else {
        diagram[pairIdentifier] = std::make_tuple(
          vertexId1, (BNodeType)nodeType1, vertexId2, (BNodeType)nodeType2,
          persistence, pairType, birth, coordsBirth[0], coordsBirth[1],
          coordsBirth[2], death, coordsDeath[0], coordsDeath[1],
          coordsDeath[2]);
      }
    }
    if(pairIdentifier >= pairingsSize) {
      nbNonCompact++;
      if(nbNonCompact == 0) {
        std::stringstream msg;
        msg << "Diagram pair identifiers "
            << "must be compact (not exceed the diagram size). " << std::endl;
        this->printWrn(msg.str());
      }
    }
  }

  if(nbNonCompact > 0) {
    {
      std::stringstream msg;
      msg << "Missed " << nbNonCompact << " pairs due to non-compactness."
          << std::endl;
      this->printWrn(msg.str());
    }
  }

  return max_dimension;
}
