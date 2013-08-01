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
#include <vtkMatrix4x4.h>

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

  void AddDavinciPrimitives(const vtkMatrix4x4& baseFrame,
                            const double* q_passive,
                            const double* q_active);

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
};

#endif
