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

// Qt includes
#include <QtPlugin>

// AutoPortPlacement Logic includes
#include <vtkSlicerAutoPortPlacementLogic.h>

// AutoPortPlacement includes
#include "qSlicerAutoPortPlacementModule.h"
#include "qSlicerAutoPortPlacementModuleWidget.h"

//-----------------------------------------------------------------------------
Q_EXPORT_PLUGIN2(qSlicerAutoPortPlacementModule, qSlicerAutoPortPlacementModule);

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerAutoPortPlacementModulePrivate
{
public:
  qSlicerAutoPortPlacementModulePrivate();
};

//-----------------------------------------------------------------------------
// qSlicerAutoPortPlacementModulePrivate methods

//-----------------------------------------------------------------------------
qSlicerAutoPortPlacementModulePrivate
::qSlicerAutoPortPlacementModulePrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerAutoPortPlacementModule methods

//-----------------------------------------------------------------------------
qSlicerAutoPortPlacementModule
::qSlicerAutoPortPlacementModule(QObject* _parent)
  : Superclass(_parent)
  , d_ptr(new qSlicerAutoPortPlacementModulePrivate)
{
}

//-----------------------------------------------------------------------------
qSlicerAutoPortPlacementModule::~qSlicerAutoPortPlacementModule()
{
}

//-----------------------------------------------------------------------------
QString qSlicerAutoPortPlacementModule::helpText()const
{
  return "This is a loadable module bundled in an extension";
}

//-----------------------------------------------------------------------------
QString qSlicerAutoPortPlacementModule::acknowledgementText()const
{
  return "This work was was partially funded by NIH grant 3P41RR013218-12S1";
}

//-----------------------------------------------------------------------------
QStringList qSlicerAutoPortPlacementModule::contributors()const
{
  QStringList moduleContributors;
  moduleContributors << QString("Jean-Christophe Fillion-Robin (Kitware)");
  return moduleContributors;
}

//-----------------------------------------------------------------------------
QIcon qSlicerAutoPortPlacementModule::icon()const
{
  return QIcon(":/Icons/AutoPortPlacement.png");
}

//-----------------------------------------------------------------------------
QStringList qSlicerAutoPortPlacementModule::categories() const
{
  return QStringList() << "Examples";
}

//-----------------------------------------------------------------------------
QStringList qSlicerAutoPortPlacementModule::dependencies() const
{
  return QStringList();
}

//-----------------------------------------------------------------------------
void qSlicerAutoPortPlacementModule::setup()
{
  this->Superclass::setup();
}

//-----------------------------------------------------------------------------
qSlicerAbstractModuleRepresentation * qSlicerAutoPortPlacementModule
::createWidgetRepresentation()
{
  return new qSlicerAutoPortPlacementModuleWidget;
}

//-----------------------------------------------------------------------------
vtkMRMLAbstractLogic* qSlicerAutoPortPlacementModule::createLogic()
{
  return vtkSlicerAutoPortPlacementLogic::New();
}
