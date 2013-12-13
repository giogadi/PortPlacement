/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

// .NAME vtkSlicerAutoPortPlacementLogic - slicer logic class for volumes manipulation
// .SECTION Description
// This class manages the logic associated with reading, saving,
// and changing propertied of the volumes


#ifndef __vtkSlicerAutoPortPlacementLogic_h
#define __vtkSlicerAutoPortPlacementLogic_h

// Slicer includes
#include "vtkSlicerModuleLogic.h"

// MRML includes

// VTK includes
#include <vtkCylinderSource.h>
#include <vtkSphereSource.h>
#include <vtkMRMLLinearTransformNode.h>

// STD includes
#include <cstdlib>

#include "vtkSlicerAutoPortPlacementModuleLogicExport.h"

// Forward declarations
class DavinciKinematics;


/// \ingroup Slicer_QtModules_ExtensionTemplate
class VTK_SLICER_AUTOPORTPLACEMENT_MODULE_LOGIC_EXPORT vtkSlicerAutoPortPlacementLogic :
  public vtkSlicerModuleLogic
{
public:

  static vtkSlicerAutoPortPlacementLogic *New();
  vtkTypeMacro(vtkSlicerAutoPortPlacementLogic, vtkSlicerModuleLogic);
  void PrintSelf(ostream& os, vtkIndent indent);

  void SetPassiveLeftJoint(unsigned jointIdx, double value);
  void SetPassiveRightJoint(unsigned jointIdx, double value);

  double GetPassiveLeftJoint(unsigned jointIdx) const;
  double GetPassiveRightJoint(unsigned jointIdx) const;

  double GetPassiveJointMin(unsigned idx) const;
  double GetPassiveJointMax(unsigned idx) const;

  void RenderRobot();

  void FindFeasiblePlan(vtkMRMLNode* taskFramesNode,
                        vtkMRMLNode* portCurvePointsNode,
                        vtkMRMLNode* robotBaseNode);

protected:
  vtkSlicerAutoPortPlacementLogic();
  virtual ~vtkSlicerAutoPortPlacementLogic();

  virtual void SetMRMLSceneInternal(vtkMRMLScene* newScene);
  /// Register MRML Node classes to Scene. Gets called automatically when the MRMLScene is attached to this logic class.
  virtual void RegisterNodes();
  virtual void UpdateFromMRMLScene();
  virtual void OnMRMLSceneNodeAdded(vtkMRMLNode* node);
  virtual void OnMRMLSceneNodeRemoved(vtkMRMLNode* node);
private:

  vtkSlicerAutoPortPlacementLogic(const vtkSlicerAutoPortPlacementLogic&); // Not implemented
  void operator=(const vtkSlicerAutoPortPlacementLogic&);               // Not implemented

  DavinciKinematics* Kinematics;

  /// For rendering the collision primitives of the davinci
  vtkSmartPointer<vtkCylinderSource> CylinderSource;
  vtkSmartPointer<vtkSphereSource> SphereSource;

  double RobotBaseX;
  double RobotBaseY;
  double RobotBaseZ;

  std::vector<double> LeftPassiveConfig;
  std::vector<double> RightPassiveConfig;
  std::vector<double> LeftActiveConfig;
  std::vector<double> RightActiveConfig;

  bool IsRobotInitialized;

  // vtkTransform BaseFrame;

  // class RobotNode
  // {
  // public:
  //   std::vector<vtkSmartPointer<vtkMRMLModelNode> > ModelNodes;
  //   std::vector<vtkSmartPointer<vtkPolyData> > PolyDataNodes;
  //   std::vector<vtkSmartPointer<vtkMRMLModelDisplayNode > DisplayNodes;
  //   std::vector<vtkSmartPointer<vtkMRMLLinearTransformNode> > TransformNodes;
  // }

  // RobotNode DavinciNode;
  std::vector<vtkSmartPointer<vtkMRMLLinearTransformNode> > RobotTransformNodes;

  void InitRobot();
  void UpdateRobot();
};

#endif
