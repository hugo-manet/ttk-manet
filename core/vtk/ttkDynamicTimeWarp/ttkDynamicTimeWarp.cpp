#include <ttkDynamicTimeWarp.h>

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkUnstructuredGrid.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <regex>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkDynamicTimeWarp);

/**
 * vTODO 7: Implement the filter constructor and destructor in the cpp file.
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
ttkDynamicTimeWarp::ttkDynamicTimeWarp() {
  this->setDebugMsgPrefix("DynamicTimeWarp");
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(2);
}

ttkDynamicTimeWarp::~ttkDynamicTimeWarp() {
}

/**
 * vTODO 8: Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkDynamicTimeWarp::FillInputPortInformation(int port,
                                                 vtkInformation *info) {
  if(port == 0)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkTable");
  else
    return 0;

  return 1;
}

/**
 * vTODO 9: Specify the data object type of each output port
 *
 * This method specifies in the port information object the data type of the
 * corresponding output objects. It is possible to either explicitly
 * specify a type by adding a vtkDataObject::DATA_TYPE_NAME() key:
 *
 *      info->Set( ttkAlgorithm::DATA_TYPE_NAME(), "vtkUnstructuredGrid" );
 *
 * or to pass a type of an input port to an output port by adding the
 * ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT() key (see below).
 *
 * Note: prior to the execution of the RequestData method the pipeline will
 * initialize empty output data objects based on this information.
 */
int ttkDynamicTimeWarp::FillOutputPortInformation(int port,
                                                  vtkInformation *info) {
  if(port == 0 || port == 1)
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  else
    return 0;

  return 1;
}

/**
 * TODO 10: Pass VTK data to the base code and convert base code output to VTK
 *
 * This method is called during the pipeline execution to update the
 * already initialized output data objects based on the given input
 * data objects and filter parameters.
 *
 * Note:
 *     1) The passed input data objects are validated based on the information
 *        provided by the FillInputPortInformation method.
 *     2) The output objects are already initialized based on the information
 *        provided by the FillOutputPortInformation method.
 */
int ttkDynamicTimeWarp::RequestData(vtkInformation *request,
                                    vtkInformationVector **inputVector,
                                    vtkInformationVector *outputVector) {
  ttk::Timer tm{};

  const auto input = vtkTable::GetData(inputVector[0]);
  std::vector<std::string> NonscalarFields;

  if(SelectFieldsWithRegexp) {
    // select all input columns whose name is matching the regexp
    ScalarFields.clear();
    NonscalarFields.clear();
    const auto n = input->GetNumberOfColumns();
    for(int jCol = 0; jCol < n; ++jCol) {
      const auto &name = input->GetColumnName(jCol);
      if(std::regex_match(name, std::regex(RegexpString))) {
        ScalarFields.emplace_back(name);
      } else {
        NonscalarFields.emplace_back(name);
      }
    }
  } else if(this->SplitMatrix && this->CopyRemainingDataOnPoints) {
    // TODO fill NonscalarFields with the absent
    // (both are sorted increasing so it's linear)
  }

  size_t nRows, nColumns;
  if(this->SplitMatrix) {
    size_t totalSize = static_cast<size_t>(input->GetNumberOfRows());
    if(ScalarFields.size() != totalSize) {
      this->printErr(
        "Distance matrix is not a square, yet you selected SplitMatrix");
      return 0;
    }
    if(this->SplitPivot >= totalSize) {
      this->printErr("SplitPivot selected higher than matrix size");
      return 0;
    }
    nRows = SplitPivot + 1;
    nColumns = totalSize - nRows;
  } else {
    nColumns = ScalarFields.size();
    nRows = static_cast<size_t>(input->GetNumberOfRows());
  }

  boost::numeric::ublas::matrix<double> distanceMatrix(nRows, nColumns);
  for(size_t iRow = 0; iRow < nRows; ++iRow) {
    for(size_t jCol = 0; jCol < nColumns; ++jCol) {
      if(this->SplitMatrix) {
        distanceMatrix(iRow, jCol)
          = input->GetColumnByName(ScalarFields[nRows + jCol].data())
              ->GetVariantValue(iRow)
              .ToDouble();
      } else {
        distanceMatrix(iRow, jCol)
          = input->GetColumnByName(ScalarFields[jCol].data())
              ->GetVariantValue(iRow)
              .ToDouble();
      }
    }
  }

  auto warpingPath
    = this->computeWarpingPath(distanceMatrix, this->DeletionCost);

  auto output_path = vtkUnstructuredGrid::SafeDownCast(
    outputVector->GetInformationObject(0)->Get(vtkDataObject::DATA_OBJECT()));
  auto output_matching = vtkUnstructuredGrid::SafeDownCast(
    outputVector->GetInformationObject(1)->Get(vtkDataObject::DATA_OBJECT()));

  vtkSmartPointer<vtkPoints> matchingPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkUnstructuredGrid> matchingCells
    = vtkSmartPointer<vtkUnstructuredGrid>::New();
  matchingCells->AllocateExact(
    warpingPath.size() + 1, 2 * (warpingPath.size() + 1));

  vtkSmartPointer<vtkPoints> pathPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkUnstructuredGrid> pathCells
    = vtkSmartPointer<vtkUnstructuredGrid>::New();

  // type is 0 for first curve (rows), 1 for second curve (cols)
  vtkNew<vtkIntArray> matchingPointCurveIndex{};
  matchingPointCurveIndex->SetNumberOfComponents(1);
  matchingPointCurveIndex->SetName("Curve_isCols");
  for(size_t iRow = 0; iRow < nRows; ++iRow) {
    matchingPoints->InsertNextPoint(iRow, 0, 0);
    matchingPointCurveIndex->InsertNextValue(0);
  }
  const size_t offsetForCols = nRows;
  for(size_t jCol = 0; jCol < nColumns; ++jCol) {
    matchingPoints->InsertNextPoint(jCol, 1, 0);
    matchingPointCurveIndex->InsertNextValue(1);
  }
  // type is 0 for intra-curve for rows, 1 for intra-curve for cols, 2 for
  // deletion, 3 for matching
  vtkNew<vtkIntArray> matchingType{};
  matchingType->SetNumberOfComponents(1);
  matchingType->SetName("Type");
  vtkNew<vtkDoubleArray> matchingDistance{};
  matchingDistance->SetNumberOfComponents(1);
  matchingDistance->SetName("Distance");

  for(size_t iRow = 1; iRow < nRows; ++iRow) {
    vtkIdType coords[2] = {iRow - 1, iRow};
    matchingCells->InsertNextCell(VTK_LINE, 2, coords);
    matchingType->InsertNextValue(0);
    if(this->SplitMatrix)
      matchingDistance->InsertNextValue(
        input->GetColumnByName(ScalarFields[iRow].data())
          ->GetVariantValue(iRow - 1)
          .ToDouble());
    else
      matchingDistance->InsertNextValue(0.);
  }
  for(size_t jCol = 1; jCol < nColumns; ++jCol) {
    vtkIdType coords[2] = {offsetForCols + jCol - 1, offsetForCols + jCol};
    matchingCells->InsertNextCell(VTK_LINE, 2, coords);
    matchingType->InsertNextValue(1);
    if(this->SplitMatrix)
      matchingDistance->InsertNextValue(
        input->GetColumnByName(ScalarFields[offsetForCols + jCol].data())
          ->GetVariantValue(offsetForCols + jCol - 1)
          .ToDouble());
    else
      matchingDistance->InsertNextValue(0.);
  }

  vtkNew<vtkDoubleArray> pathWeight{};
  pathWeight->SetNumberOfComponents(1);
  pathWeight->SetName("Weight");

  vtkNew<vtkDoubleArray> pathDistance{};
  pathDistance->SetNumberOfComponents(1);
  pathDistance->SetName("Distance");

  // Fill the matching points and the path
  vtkIdType matchingLine[2] = {0, offsetForCols};
  vtkIdType pathLine[2] = {-1, 0};
  matchingCells->InsertNextCell(VTK_LINE, 2, matchingLine);
  matchingType->InsertNextValue(3);
  matchingDistance->InsertNextValue(distanceMatrix(0, 0));
  pathPoints->InsertNextPoint(0, 0, 0);
  for(auto [dir, iRow, jCol, weightP] : warpingPath) {
    pathLine[0]++;
    pathLine[1]++;
    pathCells->InsertNextCell(VTK_LINE, 2, pathLine);
    pathPoints->InsertNextPoint(jCol, iRow, 0);
    pathDistance->InsertNextValue(distanceMatrix(iRow, jCol));
    pathWeight->InsertNextValue(weightP);
    switch(dir) {
      case Direction::DIR_SAME_COL:
        matchingLine[0] = iRow;
        matchingType->InsertNextValue(2);
        break;
      case Direction::DIR_SAME_ROW:
        matchingLine[1] = offsetForCols + jCol;
        matchingType->InsertNextValue(2);
        break;
      case Direction::DIR_BOTH:
        matchingLine[0] = iRow;
        matchingLine[1] = offsetForCols + jCol;
        matchingType->InsertNextValue(3);
        break;
    }
    matchingCells->InsertNextCell(VTK_LINE, 2, matchingLine);
    matchingDistance->InsertNextValue(distanceMatrix(iRow, jCol));
  }

  matchingCells->SetPoints(matchingPoints);
  output_matching->ShallowCopy(matchingCells);
  output_matching->GetCellData()->AddArray(matchingType);
  output_matching->GetCellData()->AddArray(matchingDistance);
  output_matching->GetPointData()->AddArray(matchingPointCurveIndex);

  pathCells->SetPoints(pathPoints);
  output_path->ShallowCopy(pathCells);
  output_path->GetCellData()->AddArray(pathWeight);
  output_path->GetCellData()->AddArray(pathDistance);

  return 1;
}
