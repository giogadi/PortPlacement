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
#include <QDebug>

// SlicerQt includes
#include "qSlicerAutoPortPlacementModuleWidget.h"
#include "ui_qSlicerAutoPortPlacementModuleWidget.h"

// Logic includes
#include "vtkSlicerAutoPortPlacementLogic.h"

// VTK includes
#include <vtkTransform.h>

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_ExtensionTemplate
class qSlicerAutoPortPlacementModuleWidgetPrivate: public Ui_qSlicerAutoPortPlacementModuleWidget
{
public:
  qSlicerAutoPortPlacementModuleWidgetPrivate();
};

//-----------------------------------------------------------------------------
// qSlicerAutoPortPlacementModuleWidgetPrivate methods

//-----------------------------------------------------------------------------
qSlicerAutoPortPlacementModuleWidgetPrivate::qSlicerAutoPortPlacementModuleWidgetPrivate()
{
}

//-----------------------------------------------------------------------------
// qSlicerAutoPortPlacementModuleWidget methods

//-----------------------------------------------------------------------------
qSlicerAutoPortPlacementModuleWidget::qSlicerAutoPortPlacementModuleWidget(QWidget* _parent)
  : Superclass( _parent )
  , d_ptr( new qSlicerAutoPortPlacementModuleWidgetPrivate )
{
}

//-----------------------------------------------------------------------------
qSlicerAutoPortPlacementModuleWidget::~qSlicerAutoPortPlacementModuleWidget()
{
}

//-----------------------------------------------------------------------------
void qSlicerAutoPortPlacementModuleWidget::setup()
{
  Q_D(qSlicerAutoPortPlacementModuleWidget);
  d->setupUi(this);  

  this->Superclass::setup();
}

void qSlicerAutoPortPlacementModuleWidget::onPushButtonClicked()
{
  if (this->logic() == 0)
    {
    return;
    }
  vtkSlicerAutoPortPlacementLogic *portLogic = 
    vtkSlicerAutoPortPlacementLogic::SafeDownCast(this->logic());
  if (portLogic)
    {
    vtkSmartPointer<vtkTransform> T1 = vtkSmartPointer<vtkTransform>::New();
    T1->Translate(-0.5, 0.0, 0.0);
    T1->RotateZ(90.0);

    vtkSmartPointer<vtkTransform> T2 = vtkSmartPointer<vtkTransform>::New();
    T2->Translate(0.5, 0.0, 0.0);
    T2->RotateZ(-90.0);

    double qpL[] = {0.454813, -1.34757, -0.124238, -1.24044, -0.659561, -0.632044};
    double qpR[] = {0.803397, 0.781243, 1.10959, 0.997729, 1.66653, -0.315385};
    double qa[] = {0, 0, 0, 0, 0, 0};

    portLogic->AddDavinciPrimitives(*T1->GetMatrix(), qpL, qa);
    portLogic->AddDavinciPrimitives(*T2->GetMatrix(), qpR, qa);
    }
}
