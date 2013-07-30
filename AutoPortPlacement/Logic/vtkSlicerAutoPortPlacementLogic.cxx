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

//----------------------------------------------------------------------------
vtkSlicerAutoPortPlacementLogic::vtkSlicerAutoPortPlacementLogic()
{
  this->Kinematics = new DavinciKinematics("/home/luis/code/PortPlacement/AutoPortPlacement/Logic/data/davinci-parameters.xml"); // obviously debug

  this->CylinderSource = vtkSmartPointer<vtkCylinderSource>::New();
  this->CylinderSource->SetHeight(1.0);
  this->CylinderSource->SetCenter(0.0,0.0,0.0);
  this->CylinderSource->SetRadius(1.0);

  this->SphereSource = vtkSmartPointer<vtkSphereSource>::New();
  this->SphereSource->SetRadius(1.0);
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

//---------------------------------------------------------------------------
void vtkSlicerAutoPortPlacementLogic::AddDavinciPrimitives()
{
  std::vector<double> q(6,0.0);
  q[0] = 0.3;
  q[1] = 1.5;
  q[2] = -0.7;
  q[3] = 0.7;
  q[4] = 1.5;
  std::vector<Collisions::Cylisphere> cylispheres;
  Collisions::Sphere sphere;
  this->Kinematics->getPassivePrimitives(Eigen::Matrix4d::Identity(),
                                         q, &cylispheres, &sphere);

  for (std::size_t i = 0; i < cylispheres.size(); ++i)
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

    // Transform cylinder into proper pose
    vtkSmartPointer<vtkTransform> t = vtkSmartPointer<vtkTransform>::New();
    
    Eigen::Vector3d p1 = cylispheres[i].p1*1000;
    Eigen::Vector3d p2 = cylispheres[i].p2*1000;

    double cylinderRad = cylispheres[i].r*1000;
    Eigen::Vector3d cylinderVec = p2 - p1;
    double cylinderHeight = cylinderVec.norm();

    // Rotate
    cylinderVec.normalize();
    Eigen::Vector3d e2; e2(0) = 0.0; e2(1) = 1.0; e2(2) = 0.0;
    Eigen::Vector3d rotationAxis = e2.cross(cylinderVec);
    double crossNorm = rotationAxis.norm();
    double angle = asin(crossNorm);
    rotationAxis /= crossNorm;


    t->Translate(p1(0), p1(1), p1(2));
    t->RotateWXYZ(vtkMath::DegreesFromRadians(angle), 
                  rotationAxis(0), rotationAxis(1), rotationAxis(2));
    t->Translate(0.0, 0.5*cylinderHeight, 0.0);
    t->Scale(cylinderRad, cylinderHeight, cylinderRad);

    vtkSmartPointer<vtkMRMLLinearTransformNode> transformNode = 
      vtkSmartPointer<vtkMRMLLinearTransformNode>::New();
    transformNode->ApplyTransformMatrix(t->GetMatrix());
    this->GetMRMLScene()->AddNode(transformNode);
   
    modelNode->SetAndObserveTransformNodeID(transformNode->GetID());    
    }

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

  vtkSmartPointer<vtkTransform> t = vtkSmartPointer<vtkTransform>::New();
  Eigen::Vector3d p = 1000*sphere.p;
  double r = 1000*sphere.r;
  t->Translate(p(0), p(1), p(2));
  t->Scale(r, r, r);
  vtkSmartPointer<vtkMRMLLinearTransformNode> transformNode = 
    vtkSmartPointer<vtkMRMLLinearTransformNode>::New();
  transformNode->ApplyTransformMatrix(t->GetMatrix());
  this->GetMRMLScene()->AddNode(transformNode);
  
  modelNode->SetAndObserveTransformNodeID(transformNode->GetID());
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

