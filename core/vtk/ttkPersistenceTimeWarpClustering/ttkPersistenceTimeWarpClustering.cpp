#include <ttkMacros.h>
#include <ttkPersistenceTimeWarpClustering.h>
#include <ttkUtils.h>
#include <vtkFieldData.h>

using namespace std;
using namespace ttk;

vtkStandardNewMacro(ttkPersistenceTimeWarpClustering)

  ttkPersistenceTimeWarpClustering::ttkPersistenceTimeWarpClustering() {
  SetNumberOfInputPorts(1);
  SetNumberOfOutputPorts(5);
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
  if(port == 2 || port == 3 || port == 4)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  else if(port == 0 || port == 1)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet");
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

  size_t nCurves = inputVector[0]->GetNumberOfInformationObjects();
  std::vector<vtkMultiBlockDataSet *> blocks(nCurves);
  std::vector<std::vector<vtkUnstructuredGrid *>> inputDiagramGrids(nCurves);

  if(nCurves == 0) {
    this->printErr("Input is empty");
    return 0;
  }

  // number of diagrams per input block
  std::vector<size_t> nDiagOfCurve(nCurves);

  for(size_t iCurve = 0; iCurve < nCurves; ++iCurve) {
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
    for(const auto &vtu : curveGrid)
      if(vtu == nullptr) {
        this->printErr("Input diagrams are not all vtkUnstructuredGrid");
        return 0;
      }
  }

  intermediateDiagramsCurves_.clear();
  max_dimension_total_ = 0;
  for(const auto &curveGrid : inputDiagramGrids) {
    intermediateDiagramsCurves_.emplace_back();
    for(const auto &vtu : curveGrid) {
      intermediateDiagramsCurves_.back().emplace_back();
      double max_dimension = this->getPersistenceDiagram(
        intermediateDiagramsCurves_.back().back(), vtu);
      if(max_dimension_total_ < max_dimension) {
        max_dimension_total_ = max_dimension;
      }
    }
  }
  final_centroid_.clear();
  final_centroid_.emplace_back(nDiagOfCurve[0]);
  inv_clustering_.assign(nCurves, 0); // everybody in one cluster

  all_matchings_.clear();
  time_warp_.clear();
  this->executeTimeWarp(intermediateDiagramsCurves_, final_centroid_[0],
                        all_matchings_, time_warp_);

  // Set outputs
  auto outputInitialDiagrams = vtkMultiBlockDataSet::SafeDownCast(
    outputVector->GetInformationObject(0)->Get(vtkDataObject::DATA_OBJECT()));
  auto outputBarycenterCurves = vtkMultiBlockDataSet::SafeDownCast(
    outputVector->GetInformationObject(1)->Get(vtkDataObject::DATA_OBJECT()));
  auto outputMatching = vtkUnstructuredGrid::SafeDownCast(
    outputVector->GetInformationObject(2)->Get(vtkDataObject::DATA_OBJECT()));
  auto outputTimeWarp = vtkUnstructuredGrid::SafeDownCast(
    outputVector->GetInformationObject(3)->Get(vtkDataObject::DATA_OBJECT()));
  auto outputSlicePoints = vtkUnstructuredGrid::SafeDownCast(
    outputVector->GetInformationObject(4)->Get(vtkDataObject::DATA_OBJECT()));

  // outputMatching->ShallowCopy(createMatchings());
  outputInitialDiagrams->ShallowCopy(createOutputClusteredDiagrams());
  outputBarycenterCurves->ShallowCopy(createOutputCentroids());
  outputTimeWarp->ShallowCopy(createOutputTimeWarp());
  outputSlicePoints->ShallowCopy(createOutputSlicePoints());

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

  // If diagram has the diagonal (we assume it is last)
  if(*pairIdentifierScalars->GetTuple(pairingsSize - 1) == -1)
    pairingsSize -= 1;

#ifndef TTK_ENABLE_KAMIKAZE
  if(pairingsSize < 1 || !vertexIdentifierScalars || !pairIdentifierScalars
     || !nodeTypeScalars || !persistenceScalars || !extremumIndexScalars
     || !points) {
    return -2;
  }
#endif // TTK_ENABLE_KAMIKAZE

  diagram.resize(pairingsSize);
  int nbNonCompact = 0;
  double max_dimension = 0;

  // skip diagonal cell (corresponding points already dealt with)
  for(int i = 0; i < pairingsSize; ++i) {

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

        diagram[0] = std::make_tuple(
          vertexId1, CriticalType::Local_minimum, vertexId2,
          CriticalType::Local_maximum, persistence, pairType, birth,
          coordsBirth[0], coordsBirth[1], coordsBirth[2], death, coordsDeath[0],
          coordsDeath[1], coordsDeath[2]);

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

vtkSmartPointer<vtkUnstructuredGrid> createDiagram(const ttk::Diagram &diagram,
                                                   double offset_x,
                                                   double offset_y,
                                                   double offset_z) {
  vtkNew<vtkUnstructuredGrid> persistenceDiagram{};
  vtkNew<vtkPoints> points{};

  vtkNew<ttkSimplexIdTypeArray> vertexIdentifierScalars{};
  vertexIdentifierScalars->SetName(ttk::VertexScalarFieldName);

  vtkNew<vtkIntArray> nodeType{};
  nodeType->SetName("CriticalType");

  vtkNew<vtkDoubleArray> persistenceScalars{};
  persistenceScalars->SetName("Persistence");

  vtkNew<vtkIntArray> idOfPair{};
  idOfPair->SetName("PairIdentifier");

  vtkNew<vtkDoubleArray> persistenceScalarsPoint{};
  persistenceScalarsPoint->SetName("Persistence");

  vtkNew<vtkIntArray> pairType{};
  pairType->SetName("PairType");

  vtkNew<vtkFloatArray> coordsScalars{};
  coordsScalars->SetNumberOfComponents(3);
  coordsScalars->SetName("Coordinates");

  // First, add diagram points to the global input diagram
  for(size_t lPair = 0; lPair < diagram.size(); ++lPair) {
    vtkIdType ids[2];
    const DiagramTuple &t = diagram[lPair];
    double x1 = std::get<6>(t) + offset_x;
    double y1 = std::get<6>(t) + offset_y;
    double z1 = offset_z;
    float coords1[3];
    coords1[0] = std::get<7>(t);
    coords1[1] = std::get<8>(t);
    coords1[2] = std::get<9>(t);

    double x2 = std::get<6>(t) + offset_x;
    double y2 = std::get<10>(t) + offset_y;
    double z2 = offset_z;

    float coords2[3];
    coords2[0] = std::get<11>(t);
    coords2[1] = std::get<12>(t);
    coords2[2] = std::get<13>(t);

    idOfPair->InsertTuple1(lPair, lPair);
    points->InsertNextPoint(x1, y1, z1);
    coordsScalars->InsertTuple3(2 * lPair, coords1[0], coords1[1], coords1[2]);
    vertexIdentifierScalars->InsertTuple1(2 * lPair, 2 * lPair);
    const ttk::CriticalType n1Type = std::get<1>(t);
    switch(n1Type) {
      case BLocalMin:
        nodeType->InsertTuple1(2 * lPair, 0);
        break;

      case BSaddle1:
        nodeType->InsertTuple1(2 * lPair, 1);
        break;

      case BSaddle2:
        nodeType->InsertTuple1(2 * lPair, 2);
        break;

      case BLocalMax:
        nodeType->InsertTuple1(2 * lPair, 3);
        break;
      default:
        nodeType->InsertTuple1(2 * lPair, 0);
    }
    points->InsertNextPoint(x2, y2, z2);
    coordsScalars->InsertTuple3(
      2 * lPair + 1, coords2[0], coords2[1], coords2[2]);
    vertexIdentifierScalars->InsertTuple1(2 * lPair + 1, 2 * lPair + 1);
    const ttk::CriticalType n2Type = std::get<3>(t);
    switch(n2Type) {
      case BLocalMin:
        nodeType->InsertTuple1(2 * lPair + 1, 0);
        break;

      case BSaddle1:
        nodeType->InsertTuple1(2 * lPair + 1, 1);
        break;

      case BSaddle2:
        nodeType->InsertTuple1(2 * lPair + 1, 2);
        break;

      case BLocalMax:
        nodeType->InsertTuple1(2 * lPair + 1, 3);
        break;
      default:
        nodeType->InsertTuple1(2 * lPair + 1, 0);
    }

    ids[0] = 2 * lPair;
    ids[1] = 2 * lPair + 1;

    persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);
    persistenceScalars->InsertTuple1(lPair, y2 - x2);
    persistenceScalarsPoint->InsertTuple1(2 * lPair, y2 - x2);
    persistenceScalarsPoint->InsertTuple1(2 * lPair + 1, y2 - x2);
    const ttk::SimplexId type = std::get<5>(t);
    switch(type) {
      case 0:
        pairType->InsertTuple1(lPair, 0);
        break;

      case 1:
        pairType->InsertTuple1(lPair, 1);
        break;

      case 2:
        pairType->InsertTuple1(lPair, 2);
        break;
      default:
        pairType->InsertTuple1(lPair, 2);
    }
  }

  // Add diagonal
  size_t newIndex = diagram.size();
  vtkIdType ids[2] = {0, 2 * newIndex - 2};
  persistenceDiagram->InsertNextCell(VTK_LINE, 2, ids);
  persistenceScalars->InsertTuple1(newIndex, -1);
  pairType->InsertTuple1(newIndex, -1);
  idOfPair->InsertTuple1(newIndex, newIndex);

  persistenceDiagram->SetPoints(points);
  persistenceDiagram->GetCellData()->AddArray(persistenceScalars);
  persistenceDiagram->GetCellData()->AddArray(pairType);
  persistenceDiagram->GetCellData()->AddArray(idOfPair);
  persistenceDiagram->GetPointData()->AddArray(nodeType);
  persistenceDiagram->GetPointData()->AddArray(coordsScalars);
  persistenceDiagram->GetPointData()->AddArray(persistenceScalarsPoint);
  persistenceDiagram->GetPointData()->AddArray(vertexIdentifierScalars);

  return persistenceDiagram;
}

vtkSmartPointer<vtkMultiBlockDataSet>
  ttkPersistenceTimeWarpClustering::createOutputCentroids() {
  this->printMsg("Creating vtk diagrams", debug::Priority::VERBOSE);

  size_t nCentroids = final_centroid_.size();
  size_t nDiags = 0;
  for(size_t jCentroid = 0; jCentroid < nCentroids; ++jCentroid)
    nDiags += final_centroid_[jCentroid].size();

  vtkNew<vtkMultiBlockDataSet> diagramCollection{};
  diagramCollection->SetNumberOfBlocks(nDiags);

  size_t iDiagFinal = 0;
  for(size_t jCentroid = 0; jCentroid < nCentroids; ++jCentroid)
    for(size_t kDiag = 0; kDiag < final_centroid_[jCentroid].size(); ++kDiag) {
      const std::vector<DiagramTuple> &diagram
        = final_centroid_[jCentroid][kDiag];

      double x = 0., y = 0., z = 0.;
      if(DisplayMethod == 1 && Spacing != 0) {
        x += 3 * (abs(Spacing) + 0.2) * max_dimension_total_ * jCentroid;
        z += kDiag * max_dimension_total_;
      }
      diagramCollection->SetBlock(
        iDiagFinal++, createDiagram(diagram, x, y, z));
    }

  return diagramCollection;
}

vtkSmartPointer<vtkMultiBlockDataSet>
  ttkPersistenceTimeWarpClustering::createOutputClusteredDiagrams() {
  this->printMsg("Creating vtk outputs", debug::Priority::VERBOSE);

  size_t nCurves = intermediateDiagramsCurves_.size();
  size_t nDiags = 0;

  std::vector<size_t> cluster_size;
  std::vector<size_t> idxInCluster(nCurves);
  for(size_t jCurve = 0; jCurve < nCurves; ++jCurve) {
    idxInCluster[jCurve] = 0;
    nDiags += intermediateDiagramsCurves_[jCurve].size();
  }
  if(Spacing > 0) {
    for(size_t jCurve = 0; jCurve < nCurves; ++jCurve) {
      size_t c = inv_clustering_[jCurve];
      if(c + 1 > cluster_size.size()) {
        cluster_size.resize(c + 1);
        cluster_size[c] = 1;
        idxInCluster[jCurve] = 0;
      } else {
        cluster_size[c]++;
        idxInCluster[jCurve] = cluster_size[c] - 1;
      }
    }
  }

  vtkNew<vtkMultiBlockDataSet> diagramCollection{};
  diagramCollection->SetNumberOfBlocks(nDiags);

  size_t iDiagFinal = 0;
  for(size_t jCurve = 0; jCurve < nCurves; ++jCurve) {
    size_t curveLength = intermediateDiagramsCurves_[jCurve].size();
    for(size_t kDiag = 0; kDiag < curveLength; ++kDiag) {

      const auto &diagram = intermediateDiagramsCurves_[jCurve][kDiag];

      double x = 0., y = 0., z = 0.;
      size_t c = inv_clustering_[jCurve];
      if(DisplayMethod == 1 && Spacing > 0) {
        double angle = 2 * 3.1415926 * idxInCluster[jCurve] / cluster_size[c];
        x += (abs(Spacing) + .2) * 3 * max_dimension_total_ * c
             + Spacing * max_dimension_total_ * cos(angle);
        y += Spacing * max_dimension_total_ * sin(angle);
        z += kDiag * max_dimension_total_;
      }
      diagramCollection->SetBlock(
        iDiagFinal++, createDiagram(diagram, x, y, z));
    }
  }

  return diagramCollection;
}

/*
vtkSmartPointer<vtkUnstructuredGrid>
  ttkPersistenceTimeWarpClustering::createMatchings() {
  this->printMsg("Creating vtk matchings", debug::Priority::VERBOSE);
  vtkNew<vtkPoints> matchingPoints{};

  vtkNew<vtkUnstructuredGrid> matchingMesh{};

  vtkNew<vtkIntArray> idOfDiagramMatchingPoint{};
  idOfDiagramMatchingPoint->SetName("DiagramID");

  vtkNew<vtkIntArray> idOfPoint{};
  idOfPoint->SetName("PointID");

  vtkNew<vtkIntArray> idOfDiagramMatching{};
  idOfDiagramMatching->SetName("DiagramID");

  vtkNew<vtkIntArray> idOfCluster{};
  idOfCluster->SetName("ClusterID");

  vtkNew<vtkDoubleArray> cost{};
  cost->SetName("Cost");

  vtkNew<vtkIntArray> pairType{};
  pairType->SetName("PairType");

  vtkNew<vtkIntArray> matchingCount{};
  matchingCount->SetName("MatchNumber");

  std::vector<int> cluster_size;
  std::vector<int> idxInCluster(intermediateDiagramsCurves_.size());

  std::vector<int> matchings_count(final_centroid_[0].size(), 0);
  std::vector<int> count_to_good;

  for(unsigned int j = 0; j < intermediateDiagramsCurves_.size(); ++j) {
    idxInCluster[j] = 0;
  }
  // RE-Invert clusters
  if(DisplayMethod == 1 && Spacing > 0) {
    for(unsigned int j = 0; j < intermediateDiagramsCurves_.size(); ++j) {
      unsigned int c = inv_clustering_[j];
      if(c + 1 > cluster_size.size()) {
        cluster_size.resize(c + 1);
        cluster_size[c] = 1;
        idxInCluster[j] = 0;
      } else {
        cluster_size[c]++;
        idxInCluster[j] = cluster_size[c] - 1;
      }
    }
  }
  int count = 0;
  for(unsigned int j = 0; j < intermediateDiagramsCurves_.size(); ++j) {
    int c = inv_clustering_[j];
    const auto &diagram = intermediateDiagramsCurves_[j];
    std::vector<matchingType> matchings_j
      = all_matchings_[inv_clustering_[j]][j];
    for(unsigned int i = 0; i < matchings_j.size(); ++i) {

      vtkIdType ids[2];
      ids[0] = 2 * count;
      ids[1] = 2 * count + 1;
      matchingType m = matchings_j[i];
      const size_t bidder_id = std::get<0>(m);
      const size_t good_id = std::get<1>(m);

      // avoid out-of-bound accesses
      if(good_id >= matchings_count.size() || bidder_id >= diagram.size()) {
        continue;
      }

      if(NumberOfClusters == 1) {
        matchings_count[good_id] += 1;
        count_to_good.push_back(good_id);
      }

      DiagramTuple t1 = final_centroid_[c][good_id];
      double x1 = std::get<6>(t1);
      double y1 = std::get<10>(t1);
      double z1 = 0;

      DiagramTuple t2 = diagram[bidder_id];
      double x2 = std::get<6>(t2);
      double y2 = std::get<10>(t2);
      double z2 = 0; // Change 1 to j if you want to isolate the diagrams

      if(DisplayMethod == 1 && Spacing > 0) {
        double angle
          = 2 * 3.1415926 * (double)(idxInCluster[j]) / cluster_size[c];
        x1 += (abs(Spacing) + .2) * 3 * max_dimension_total_ * c;
        x2 += (abs(Spacing) + .2) * 3 * max_dimension_total_ * c
              + Spacing * max_dimension_total_ * cos(angle);
        y2 += Spacing * max_dimension_total_ * sin(angle);
      } else if(DisplayMethod == 2) {
        z2 = Spacing;
        if(intermediateDiagramsCurves_.size() == 2 and j == 0) {
          z2 = -Spacing;
        }
      }

      matchingPoints->InsertNextPoint(x1, y1, z1);
      matchingPoints->InsertNextPoint(x2, y2, z2);
      matchingMesh->InsertNextCell(VTK_LINE, 2, ids);
      idOfDiagramMatching->InsertTuple1(count, j);
      idOfCluster->InsertTuple1(count, inv_clustering_[j]);
      cost->InsertTuple1(count, std::get<2>(m));
      idOfDiagramMatchingPoint->InsertTuple1(2 * count, j);
      idOfDiagramMatchingPoint->InsertTuple1(2 * count + 1, j);
      idOfPoint->InsertTuple1(2 * count, good_id);
      idOfPoint->InsertTuple1(2 * count + 1, bidder_id);

      const ttk::SimplexId type = std::get<5>(t2);
      switch(type) {
        case 0:
          pairType->InsertTuple1(count, 0);
          break;

        case 1:
          pairType->InsertTuple1(count, 1);
          break;

        case 2:
          pairType->InsertTuple1(count, 2);
          break;
        default:
          pairType->InsertTuple1(count, 0);
      }
      count++;
    }
  }

  if(NumberOfClusters == 1 and intermediateDiagramsCurves_.size() == 2) {
    for(int i = 0; i < count; i++) {
      matchingCount->InsertTuple1(i, matchings_count[count_to_good[i]]);
    }
  }

  matchingMesh->SetPoints(matchingPoints);
  matchingMesh->GetPointData()->AddArray(idOfDiagramMatchingPoint);
  matchingMesh->GetPointData()->AddArray(idOfPoint);
  matchingMesh->GetCellData()->AddArray(idOfDiagramMatching);
  matchingMesh->GetCellData()->AddArray(idOfCluster);
  matchingMesh->GetCellData()->AddArray(pairType);
  matchingMesh->GetCellData()->AddArray(cost);
  if(NumberOfClusters == 1 and intermediateDiagramsCurves_.size() == 2) {
    matchingMesh->GetCellData()->AddArray(matchingCount);
  }

  return matchingMesh;
}
// */

vtkSmartPointer<vtkUnstructuredGrid>
  ttkPersistenceTimeWarpClustering::createOutputTimeWarp() {
  this->printMsg("Creating vtk time warp", debug::Priority::VERBOSE);
  vtkNew<vtkPoints> points{};

  vtkNew<vtkUnstructuredGrid> timeWarpResult{};

  vtkNew<vtkIntArray> weightCells{};
  weightCells->SetName("Weight");

  vtkNew<vtkIntArray> sliceIdCells{};
  sliceIdCells->SetName("SliceID");

  vtkNew<vtkIntArray> idOfDiagramPoint{};
  idOfDiagramPoint->SetName("DiagramID");

  vtkNew<vtkIntArray> idOfCurvePoint{};
  idOfCurvePoint->SetName("CurveID");

  vtkNew<vtkIntArray> idOfClusterPoint{};
  idOfClusterPoint->SetName("ClusterID");

  std::vector<size_t> cluster_size;
  std::vector<size_t> idxInCluster(intermediateDiagramsCurves_.size());
  for(size_t jCurve = 0; jCurve < intermediateDiagramsCurves_.size();
      ++jCurve) {
    idxInCluster[jCurve] = 0;
  }

  if(Spacing > 0) {
    for(size_t jCurve = 0; jCurve < intermediateDiagramsCurves_.size();
        ++jCurve) {
      size_t c = inv_clustering_[jCurve];
      if(c + 1 > cluster_size.size()) {
        cluster_size.resize(c + 1);
        cluster_size[c] = 1;
        idxInCluster[jCurve] = 0;
      } else {
        cluster_size[c]++;
        idxInCluster[jCurve] = cluster_size[c] - 1;
      }
    }
  }

  size_t count = 0;
  // Compute curve point indexes and build them
  std::vector<size_t> startOfCurve = {0};
  for(size_t jCurve = 0; jCurve < intermediateDiagramsCurves_.size();
      ++jCurve) {
    for(size_t kDiag = 0; kDiag < intermediateDiagramsCurves_[jCurve].size();
        ++kDiag) {
      double x = 0;
      double y = 0;
      if(DisplayMethod == 1 && Spacing > 0) {
        double angle = 2 * 3.1415926 * (double)(idxInCluster[jCurve])
                       / cluster_size[inv_clustering_[jCurve]];
        x += (abs(Spacing) + .2) * 3 * max_dimension_total_
               * inv_clustering_[jCurve]
             + Spacing * max_dimension_total_ * cos(angle);
        y += Spacing * max_dimension_total_ * sin(angle);
      }
      double z = kDiag * max_dimension_total_;
      points->InsertNextPoint(x, y, z);
      idOfDiagramPoint->InsertNextValue(count++);
      idOfCurvePoint->InsertNextValue(jCurve);
      idOfClusterPoint->InsertNextValue(inv_clustering_[jCurve]);
    }
    startOfCurve.push_back(startOfCurve.back()
                           + intermediateDiagramsCurves_[jCurve].size());
  }
  // Last isn't a curve, it's a centroid
  // Compute centroid points indexes and build them
  std::vector<size_t> startOfCentroid = {startOfCurve.back()};
  startOfCurve.pop_back();
  for(size_t iCentroid = 0; iCentroid < final_centroid_.size(); ++iCentroid) {
    for(size_t kDiag = 0; kDiag < final_centroid_[iCentroid].size(); ++kDiag) {
      double x = 0;
      double y = 0;
      if(DisplayMethod == 1 && Spacing > 0)
        x += (abs(Spacing) + .2) * 3 * max_dimension_total_ * iCentroid;
      double z = kDiag * max_dimension_total_;
      points->InsertNextPoint(x, y, z);
      idOfDiagramPoint->InsertNextValue(count++);
      idOfCurvePoint->InsertNextValue(intermediateDiagramsCurves_.size()
                                      + iCentroid);
      idOfClusterPoint->InsertNextValue(iCentroid);
    }
    startOfCentroid.push_back(startOfCentroid.back()
                              + final_centroid_[iCentroid].size());
  }
  startOfCentroid.pop_back();

  for(size_t iCentroid = 0; iCentroid < time_warp_.size(); ++iCentroid)
    for(size_t jCurve = 0; jCurve < time_warp_[iCentroid].size(); ++jCurve)
      for(auto &[kCentroidID, lCurveID, weight] :
          time_warp_[iCentroid][jCurve]) {
        vtkIdType ids[2] = {kCentroidID + startOfCentroid[iCentroid],
                            lCurveID + startOfCurve[jCurve]};
        timeWarpResult->InsertNextCell(VTK_LINE, 2, ids);
        weightCells->InsertNextValue(weight);
        sliceIdCells->InsertNextValue(kCentroidID);
      }

  timeWarpResult->SetPoints(points);
  timeWarpResult->GetPointData()->AddArray(idOfDiagramPoint);
  timeWarpResult->GetPointData()->AddArray(idOfCurvePoint);
  timeWarpResult->GetPointData()->AddArray(idOfClusterPoint);
  timeWarpResult->GetCellData()->AddArray(weightCells);
  timeWarpResult->GetCellData()->AddArray(sliceIdCells);

  return timeWarpResult;
}

vtkSmartPointer<vtkUnstructuredGrid>
  ttkPersistenceTimeWarpClustering::createOutputSlicePoints() {
  this->printMsg("Creating vtk slice points", debug::Priority::VERBOSE);
  vtkNew<vtkPoints> points{};

  vtkNew<vtkUnstructuredGrid> slicePointsResult{};

  vtkNew<vtkIntArray> idOfDiagramPoint{};
  idOfDiagramPoint->SetName("SliceID");

  vtkNew<vtkIntArray> idOfClusterPoint{};
  idOfClusterPoint->SetName("ClusterID");

  std::vector<size_t> offsetForCurve = {0};
  for(auto &curve : intermediateDiagramsCurves_)
    offsetForCurve.push_back(offsetForCurve.back() + curve.size());

  for(size_t iCentroid = 0; iCentroid < time_warp_.size(); ++iCentroid) {
    const size_t nbSlices = final_centroid_[iCentroid].size();
    std::vector<std::vector<std::tuple<size_t, size_t>>> sliceMembers(nbSlices);
    for(size_t jCurve = 0; jCurve < time_warp_[iCentroid].size(); ++jCurve)
      for(auto &[kCentroidID, lCurveID, w] : time_warp_[iCentroid][jCurve])
        sliceMembers[kCentroidID].push_back({jCurve, lCurveID});
    for(size_t kSlice = 0; kSlice < nbSlices; ++kSlice) {
      for(auto &[j1, l1] : sliceMembers[kSlice]) {
        points->InsertNextPoint(
          offsetForCurve[j1] + l1, offsetForCurve.back() + kSlice, 0);
        idOfDiagramPoint->InsertNextValue(kSlice);
        idOfClusterPoint->InsertNextValue(iCentroid);
        for(auto &[j2, l2] : sliceMembers[kSlice]) {
          points->InsertNextPoint(
            offsetForCurve[j1] + l1, offsetForCurve[j2] + l2, 0);
          idOfDiagramPoint->InsertNextValue(kSlice);
          idOfClusterPoint->InsertNextValue(iCentroid);
        }
      }
    }
  }
  slicePointsResult->SetPoints(points);
  slicePointsResult->GetPointData()->AddArray(idOfDiagramPoint);
  slicePointsResult->GetPointData()->AddArray(idOfClusterPoint);

  return slicePointsResult;
}
