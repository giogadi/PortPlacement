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

// AutoPortPlacement Logic includes
#include "vtkSlicerAutoPortPlacementLogic.h"

// Kinematics includes
#include <davinci-kinematics/davinci.h>

// MRML includes
#include <vtkMRMLModelNode.h>
#include <vtkMRMLModelDisplayNode.h>
#include <vtkMRMLLinearTransformNode.h>

// VTK includes
#include <vtkNew.h>
#include <vtkTransform.h>
#include <vtkMath.h>

// STD includes
#include <cassert>

//----------------------------------------------------------------------------
vtkStandardNewMacro(vtkSlicerAutoPortPlacementLogic);

namespace
{
  // inline vtkTransToEigenMatrix(const vtkTransform& matrixVTK,
  //                              Eigen::Matrix4d* matrixEigen)
  // {
  //   for (unsigned i = 0; i < 4; ++i)
  //     for (unsigned j = 0; j < 4; ++j)
  //       (*matrixEigen)(i,j) = (*(matrixVTK.GetMatrix()))[i][j];
  // }

  void cylisphereToTransformNode(const Collisions::Cylisphere& c,
                                 vtkSmartPointer<vtkMRMLLinearTransformNode> tNode)
  {
    Eigen::Vector3d p1 = c.p1*1000;
    Eigen::Vector3d p2 = c.p2*1000;

    double cylinderRad = c.r*1000;
    Eigen::Vector3d cylinderVec = p2 - p1;
    double cylinderHeight = cylinderVec.norm();

    // Rotate
    cylinderVec.normalize();
    Eigen::Vector3d e2; e2(0) = 0.0; e2(1) = 1.0; e2(2) = 0.0;
    Eigen::Vector3d rotationAxis = e2.cross(cylinderVec);
    double angle = acos(e2.dot(cylinderVec));
    rotationAxis.normalize();

    vtkSmartPointer<vtkMatrix4x4> m = vtkSmartPointer<vtkMatrix4x4>::New();
    tNode->GetMatrixTransformToWorld(m.GetPointer());
    vtkSmartPointer<vtkTransform> t = vtkSmartPointer<vtkTransform>::New();
    t->SetMatrix(m);
    t->Inverse();
    tNode->ApplyTransformMatrix(t->GetMatrix());

    t->Identity();
    t->Translate(p1(0), p1(1), p1(2));
    t->RotateWXYZ(vtkMath::DegreesFromRadians(angle), 
                  rotationAxis(0), rotationAxis(1), rotationAxis(2));
    t->Translate(0.0, 0.5*cylinderHeight, 0.0);
    t->Scale(cylinderRad, cylinderHeight, cylinderRad);

    tNode->ApplyTransformMatrix(t->GetMatrix());
  }

  void sphereToTransformNode(const Collisions::Sphere& s,
                             vtkSmartPointer<vtkMRMLLinearTransformNode> tNode)
  {
    vtkSmartPointer<vtkMatrix4x4> m = vtkSmartPointer<vtkMatrix4x4>::New();
    tNode->GetMatrixTransformToWorld(m.GetPointer());
    vtkSmartPointer<vtkTransform> t = vtkSmartPointer<vtkTransform>::New();
    t->SetMatrix(m);
    t->Inverse();
    tNode->ApplyTransformMatrix(t->GetMatrix());

    t->Identity();
    Eigen::Vector3d p = 1000*s.p;
    double r = 1000*s.r;
    t->Translate(p(0), p(1), p(2));
    t->Scale(r, r, r);

    tNode->ApplyTransformMatrix(t->GetMatrix());
  }
}

//----------------------------------------------------------------------------
vtkSlicerAutoPortPlacementLogic::vtkSlicerAutoPortPlacementLogic() :
  RobotBaseX(0.0),
  RobotBaseY(0.0),
  LeftPassiveConfig(6),
  RightPassiveConfig(6),
  LeftActiveConfig(6),
  RightActiveConfig(6),
  IsRobotInitialized(false)
{
  this->Kinematics = new DavinciKinematics();

  this->CylinderSource = vtkSmartPointer<vtkCylinderSource>::New();
  this->CylinderSource->SetHeight(1.0);
  this->CylinderSource->SetCenter(0.0,0.0,0.0);
  this->CylinderSource->SetRadius(1.0);

  this->SphereSource = vtkSmartPointer<vtkSphereSource>::New();
  this->SphereSource->SetRadius(1.0);

  this->Kinematics->getDefaultPassiveConfig(&(this->LeftPassiveConfig));
  this->Kinematics->getDefaultPassiveConfig(&(this->RightPassiveConfig));
  this->Kinematics->getDefaultActiveConfig(&(this->LeftActiveConfig));
  this->Kinematics->getDefaultActiveConfig(&(this->RightActiveConfig));
}

//----------------------------------------------------------------------------
vtkSlicerAutoPortPlacementLogic::~vtkSlicerAutoPortPlacementLogic()
{
  delete this->Kinematics;
}

//----------------------------------------------------------------------------
void vtkSlicerAutoPortPlacementLogic::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//----------------------------------------------------------------------------
void vtkSlicerAutoPortPlacementLogic::RenderRobot()
{
  if (this->IsRobotInitialized)
    this->UpdateRobot();
  else
    this->InitRobot();
}

//----------------------------------------------------------------------------
void vtkSlicerAutoPortPlacementLogic::InitRobot()
{
  std::vector<Collisions::Cylisphere> cylispheres;
  std::vector<Collisions::Sphere> spheres(2);
  Eigen::Matrix4d baseFrameL = Eigen::Matrix4d::Identity();
  baseFrameL(0,0) = -1.0;
  baseFrameL(1,1) = -1.0;
  baseFrameL(0,3) = this->RobotBaseX;
  baseFrameL(1,3) = this->RobotBaseY;
  this->Kinematics->getPassivePrimitives(baseFrameL,
                                         this->LeftPassiveConfig,
                                         &cylispheres,
                                         &(spheres[0]));
  Eigen::Matrix4d baseFrameR = Eigen::Matrix4d::Identity();
  baseFrameR(0,3) = this->RobotBaseX;
  baseFrameR(1,3) = this->RobotBaseY;
  this->Kinematics->getPassivePrimitives(baseFrameR,
                                         this->RightPassiveConfig,
                                         &cylispheres,
                                         &(spheres[1]));
  this->Kinematics->getExtraCylispheres(this->Kinematics->passiveFK(baseFrameL,
                                                                    this->LeftPassiveConfig),
                                        this->LeftActiveConfig,
                                        &cylispheres);
  this->Kinematics->getExtraCylispheres(this->Kinematics->passiveFK(baseFrameR,
                                                                    this->RightPassiveConfig),
                                        this->RightActiveConfig,
                                        &cylispheres);

  for (unsigned i = 0; i < cylispheres.size(); ++i)
    {
    vtkSmartPointer<vtkMRMLModelNode> modelNode = vtkSmartPointer<vtkMRMLModelNode>::New();
    modelNode->SetName(this->GetMRMLScene()->GenerateUniqueName("Link").c_str());
    
    vtkSmartPointer<vtkPolyData> polyData = this->CylinderSource->GetOutput();
    modelNode->SetAndObservePolyData(polyData);

    this->GetMRMLScene()->AddNode(modelNode);
    
    vtkSmartPointer<vtkMRMLModelDisplayNode> displayNode = 
      vtkSmartPointer<vtkMRMLModelDisplayNode>::New();
    displayNode->SetColor(0.0,1.0,1.0); // cyan
    this->GetMRMLScene()->AddNode(displayNode);
    modelNode->SetAndObserveDisplayNodeID(displayNode->GetID());

    vtkSmartPointer<vtkMRMLLinearTransformNode> transformNode = 
      vtkSmartPointer<vtkMRMLLinearTransformNode>::New();
    cylisphereToTransformNode(cylispheres[i],
                              transformNode);
    this->GetMRMLScene()->AddNode(transformNode);   
    modelNode->SetAndObserveTransformNodeID(transformNode->GetID());

    this->RobotTransformNodes.push_back(transformNode);
    }

  for (unsigned i = 0; i < spheres.size(); ++i)
    {
    vtkSmartPointer<vtkMRMLModelNode> modelNode = vtkSmartPointer<vtkMRMLModelNode>::New();
    modelNode->SetName(this->GetMRMLScene()->GenerateUniqueName("Link").c_str());
    vtkSmartPointer<vtkPolyData> polyData = this->SphereSource->GetOutput();
    modelNode->SetAndObservePolyData(polyData);

    this->GetMRMLScene()->AddNode(modelNode);
    
    vtkSmartPointer<vtkMRMLModelDisplayNode> displayNode = 
      vtkSmartPointer<vtkMRMLModelDisplayNode>::New();
    displayNode->SetColor(0.0,1.0,1.0); // cyan
    this->GetMRMLScene()->AddNode(displayNode);
    modelNode->SetAndObserveDisplayNodeID(displayNode->GetID());

    vtkSmartPointer<vtkMRMLLinearTransformNode> transformNode = 
      vtkSmartPointer<vtkMRMLLinearTransformNode>::New();
    sphereToTransformNode(spheres[i], transformNode);
    this->GetMRMLScene()->AddNode(transformNode);  
    modelNode->SetAndObserveTransformNodeID(transformNode->GetID());

    this->RobotTransformNodes.push_back(transformNode);
    }

  this->IsRobotInitialized = true;
}

//---------------------------------------------------------------------------
void vtkSlicerAutoPortPlacementLogic::UpdateRobot()
{
  std::vector<Collisions::Cylisphere> cylispheres;
  std::vector<Collisions::Sphere> spheres(2);
  Eigen::Matrix4d baseFrameL = Eigen::Matrix4d::Identity();
  baseFrameL(0,0) = -1.0;
  baseFrameL(1,1) = -1.0;
  baseFrameL(0,3) = this->RobotBaseX;
  baseFrameL(1,3) = this->RobotBaseY;
  this->Kinematics->getPassivePrimitives(baseFrameL,
                                         this->LeftPassiveConfig,
                                         &cylispheres,
                                         &(spheres[0]));
  Eigen::Matrix4d baseFrameR = Eigen::Matrix4d::Identity();
  baseFrameR(0,3) = this->RobotBaseX;
  baseFrameR(1,3) = this->RobotBaseY;
  this->Kinematics->getPassivePrimitives(baseFrameR,
                                         this->RightPassiveConfig,
                                         &cylispheres,
                                         &(spheres[1]));
  this->Kinematics->getExtraCylispheres(this->Kinematics->passiveFK(baseFrameL,
                                                                    this->LeftPassiveConfig),
                                        this->LeftActiveConfig,
                                        &cylispheres);
  this->Kinematics->getExtraCylispheres(this->Kinematics->passiveFK(baseFrameR,
                                                                    this->RightPassiveConfig),
                                        this->RightActiveConfig,
                                        &cylispheres);

  for (unsigned i = 0; i < cylispheres.size(); ++i)
    {
    cylisphereToTransformNode(cylispheres[i],
                              this->RobotTransformNodes[i]);
    }
  for (unsigned i = 0; i < spheres.size(); ++i)
    {
    sphereToTransformNode(spheres[i],
                          this->RobotTransformNodes[cylispheres.size()+i]);
    }
}

//---------------------------------------------------------------------------
void vtkSlicerAutoPortPlacementLogic::SetPassiveLeftJoint(unsigned jointIdx, double value)
{
  this->LeftPassiveConfig[jointIdx] = value;
}

//---------------------------------------------------------------------------
void vtkSlicerAutoPortPlacementLogic::SetPassiveRightJoint(unsigned jointIdx, double value)
{
  this->RightPassiveConfig[jointIdx] = value;
}

//---------------------------------------------------------------------------
double vtkSlicerAutoPortPlacementLogic::GetPassiveLeftJoint(unsigned jointIdx) const
{
  return this->LeftPassiveConfig[jointIdx];
}

//---------------------------------------------------------------------------
double vtkSlicerAutoPortPlacementLogic::GetPassiveRightJoint(unsigned jointIdx) const
{
  return this->RightPassiveConfig[jointIdx];
}

//---------------------------------------------------------------------------
double vtkSlicerAutoPortPlacementLogic::GetPassiveJointMin(unsigned idx) const
{
  return this->Kinematics->getPassiveJointMin(idx);
}

//---------------------------------------------------------------------------
double vtkSlicerAutoPortPlacementLogic::GetPassiveJointMax(unsigned idx) const
{
  return this->Kinematics->getPassiveJointMax(idx);
}

//---------------------------------------------------------------------------
void vtkSlicerAutoPortPlacementLogic::SetMRMLSceneInternal(vtkMRMLScene * newScene)
{
  vtkNew<vtkIntArray> events;
  events->InsertNextValue(vtkMRMLScene::NodeAddedEvent);
  events->InsertNextValue(vtkMRMLScene::NodeRemovedEvent);
  events->InsertNextValue(vtkMRMLScene::EndBatchProcessEvent);
  this->SetAndObserveMRMLSceneEventsInternal(newScene, events.GetPointer());
}

//-----------------------------------------------------------------------------
void vtkSlicerAutoPortPlacementLogic::RegisterNodes()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerAutoPortPlacementLogic::UpdateFromMRMLScene()
{
  assert(this->GetMRMLScene() != 0);
}

//---------------------------------------------------------------------------
void vtkSlicerAutoPortPlacementLogic
::OnMRMLSceneNodeAdded(vtkMRMLNode* vtkNotUsed(node))
{
}

//---------------------------------------------------------------------------
void vtkSlicerAutoPortPlacementLogic
::OnMRMLSceneNodeRemoved(vtkMRMLNode* vtkNotUsed(node))
{
}

