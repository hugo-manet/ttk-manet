/// \ingroup vtk
/// \class ttkPointDataSelector
/// \author Guillaume Favelier <guillaume.favelier@lip6.fr>
/// \date December 2017
///
/// \brief TTK VTK-filter that selects scalar fields on input with shallow copy.
///
/// \param Input Input scalar field (vtkDataSet)
/// \param Output Output scalar field (vtkDataSet)
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
#pragma once

#include <array>
#include <limits>
#include <string>

// VTK includes
#include <vtkCharArray.h>
#include <vtkDataArray.h>
#include <vtkDataArraySelection.h>
#include <vtkDataSet.h>
#include <vtkDataSetAlgorithm.h>
#include <vtkDoubleArray.h>
#include <vtkFiltersCoreModule.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkInformation.h>
#include <vtkIntArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedShortArray.h>

// VTK Module
#include <ttkPointDataSelectorModule.h>

// ttk code includes
#include <ttkAlgorithm.h>

class TTKPOINTDATASELECTOR_EXPORT ttkPointDataSelector : public ttkAlgorithm {

public:
  static ttkPointDataSelector *New();
  vtkTypeMacro(ttkPointDataSelector, ttkAlgorithm);

  vtkSetMacro(RegexpString, std::string);

  void SetRangeId(int data0, int data1) {
    RangeId[0] = data0;
    RangeId[1] = data1;
    Modified();
  }
  int *GetRangeId() {
    return RangeId.data();
  }

  vtkSetMacro(RenameSelected, bool);
  vtkGetMacro(RenameSelected, bool);

  vtkSetMacro(SelectedFieldName, std::string);
  vtkGetMacro(SelectedFieldName, std::string);

  void AddScalarField(std::string s) {
    SelectedFields.emplace_back(s);
    Modified();
  }

  void ClearScalarFields() {
    SelectedFields.clear();
    Modified();
  }

  vtkDataArraySelection *GetRangeIds() {
    vtkDataArraySelection *arr = vtkDataArraySelection::New();
    arr->SetArraySetting("0", true);
    arr->SetArraySetting(
      std::to_string(AvailableFields.size() - 1).c_str(), true);
    return arr;
  }

protected:
  ttkPointDataSelector();

  int RequestInformation(vtkInformation *request,
                         vtkInformationVector **inputVector,
                         vtkInformationVector *outputVector) override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  void FillAvailableFields(vtkDataSet *input);

private:
  bool RenameSelected{false};
  std::string SelectedFieldName{"SelectedField"};
  std::vector<std::string> SelectedFields{};
  std::vector<std::string> AvailableFields{};
  std::string RegexpString{".*"};
  vtkSmartPointer<vtkDataArray> localFieldCopy_{};
  std::array<int, 2> RangeId{0, std::numeric_limits<int>::max()};

  int doIt(vtkDataSet *input, vtkDataSet *output);
};
