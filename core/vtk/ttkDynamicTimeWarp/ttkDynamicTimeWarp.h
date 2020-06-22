/// vTODO 4: Provide your information
///
/// \ingroup vtk
/// \class ttkDynamicTimeWarp
/// \author Hugo Manet <hugo.manet@ens.fr>
/// \date 2020-06-18
///
/// \brief TTK VTK-filter that wraps the ttk::DynamicTimeWarp module.
/// TODO
/// This VTK filter uses the ttk::DynamicTimeWarp module to compute the bounding box
/// of a vtkDataSet, which is returned as a vtkUnstructuredGrid.
///
/// \param Input vtkDataSet whose bounding box will be computed.
/// \param Output vtkUnstructuredGrid that corresponds to bounding box of the
/// input.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::DynamicTimeWarp
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkDynamicTimeWarpModule.h>

// VTK Includes
#include <ttkAlgorithm.h>

// TTK Base Includes
#include <DynamicTimeWarp.h>

class TTKDYNAMICTIMEWARP_EXPORT ttkDynamicTimeWarp
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::DynamicTimeWarp // and we inherit from the base class
{
private:
  /**
   * vTODO 5: Add all filter parameters only as private member variables and
   *         initialize them here.
   */
  std::string OutputPathName{"WarpingPath"};
  std::string OutputMatchingName{"TimeWarpMatching"};
  /** Configures the weights on the path 
   */
  double DeletionCost{0.5};

  /** Copy filter mechanism from MatrixToHeatMap
   */
  bool SelectFieldsWithRegexp{false};
  std::string RegexpString{".*"};
  std::vector<std::string> ScalarFields;

public:
  /**
   * vTODO 6: Automatically generate getters and setters of filter
   *         parameters via vtkMacros.
   */
  vtkSetMacro(OutputPathName, std::string);
  vtkGetMacro(OutputPathName, std::string);
  vtkSetMacro(OutputMatchingName, std::string);
  vtkGetMacro(OutputMatchingName, std::string);
  vtkSetMacro(DeletionCost, double);
  vtkGetMacro(DeletionCost, double);

  /** TODO maybe put these two in .cpp ?
   */
  void SetScalarFields(const std::string &s) {
    ScalarFields.emplace_back(s);
    Modified();
  }
  void ClearScalarFields() {
    ScalarFields.clear();
    Modified();
  }

  vtkSetMacro(SelectFieldsWithRegexp, bool);
  vtkGetMacro(SelectFieldsWithRegexp, bool);

  vtkSetMacro(RegexpString, std::string);
  vtkGetMacro(RegexpString, std::string);

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkDynamicTimeWarp *New();
  vtkTypeMacro(ttkDynamicTimeWarp, ttkAlgorithm);

protected:
  /**
   * TODO 7: Implement the filter constructor and destructor
   *         (see cpp file)
   */
  ttkDynamicTimeWarp();
  ~ttkDynamicTimeWarp() override;

  /**
   * TODO 8: Specify the input data type of each input port
   *         (see cpp file)
   */
  int FillInputPortInformation(int port, vtkInformation *info) override;

  /**
   * TODO 9: Specify the data object type of each output port
   *         (see cpp file)
   */
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  /**
   * TODO 10: Pass VTK data to the base code and convert base code output to VTK
   *          (see cpp file)
   */
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
};
