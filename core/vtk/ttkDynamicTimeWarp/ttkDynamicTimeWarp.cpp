#include <ttkDynamicTimeWarp.h>

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
int ttkDynamicTimeWarp::FillInputPortInformation(int port, vtkInformation *info) {
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
int ttkDynamicTimeWarp::FillOutputPortInformation(int port, vtkInformation *info) {
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

  if(SelectFieldsWithRegexp) {
    // select all input columns whose name is matching the regexp
    ScalarFields.clear();
    const auto n = input->GetNumberOfColumns();
    for(int jCol = 0; jCol < n; ++jCol) {
      const auto &name = input->GetColumnName(jCol);
      if(std::regex_match(name, std::regex(RegexpString))) {
        ScalarFields.emplace_back(name);
      }
    }
  }

  const auto nColumns = ScalarFields.size();
  const auto nRows = static_cast<size_t>(input->GetNumberOfRows());

  boost::numeric::ublas::matrix<double> distanceMatrix(nRows, nColumns);
  for(size_t iRow = 0; iRow < nRows; ++iRow) {
    for(size_t jCol = 0; jCol < nColumns; ++jCol) {
      distanceMatrix(iRow, jCol)
        = input->GetColumnByName(ScalarFields[jCol].data())
            ->GetVariantValue(iRow)
            .ToDouble();
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

  vtkSmartPointer<vtkPoints> pathPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkUnstructuredGrid> pathCells
    = vtkSmartPointer<vtkUnstructuredGrid>::New();

  // Fill the matching points and the path
  // TODO update the cells, and their infos
  int iRow = 0, jCol = 0;
  matchingPoints->InsertNextPoint(iRow, 0, 0); // y=0 : curve along rows
  matchingPoints->InsertNextPoint(jCol, 1, 0); // y=1 : curve along cols
  pathPoints->InsertNextPoint(iRow, jCol, 0);
  int kLastRow = 0, kLastCol = 1, kLastPoint = 1;
  for(auto dir : warpingPath) {
    switch(dir) {
      case Direction::DIR_SAME_COL:
        matchingPoints->InsertNextPoint(++iRow, 0, 0);
        pathPoints->InsertNextPoint(iRow, jCol, 0);
        // TODO connect two last path points with algo weight
        // TODO connect kLastRow and kPoint along same curve
        // TODO connect kLastCol and kPoint with a matching
        kLastRow = ++kLastPoint; // new point is a row
        break;
      case Direction::DIR_SAME_ROW:
        matchingPoints->InsertNextPoint(++jCol, 1, 0);
        pathPoints->InsertNextPoint(iRow, jCol, 0);
        // TODO connect two last path points with algo weight
        // TODO connect kLastCol and kPoint along same curve
        // TODO connect kLastRow and kPoint with a matching
        kLastCol = ++kLastPoint; // new point is a col
        break;
      case Direction::DIR_BOTH:
        matchingPoints->InsertNextPoint(++iRow, 0, 0);
        matchingPoints->InsertNextPoint(++jCol, 1, 0);
        pathPoints->InsertNextPoint(iRow, jCol, 0);
        // TODO connect two last path points with algo weight
        // TODO connect both points along same curve
        // TODO connect both points with a matching
        kLastRow = ++kLastPoint; // first new point is a row
        kLastCol = ++kLastPoint; // second new point is a col
        break;
    }
  }

  matchingCells->SetPoints(matchingPoints);
  output_matching->ShallowCopy(matchingCells);

  pathCells->SetPoints(pathPoints);
  output_path->ShallowCopy(pathCells);

  return 1;
}
