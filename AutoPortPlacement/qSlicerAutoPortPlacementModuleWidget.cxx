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
    T1->RotateZ(180.0);

    vtkSmartPointer<vtkTransform> T2 = vtkSmartPointer<vtkTransform>::New();

    // double qpL[] = {0.5, -0.17, -0.17, -0.17, -0.17, -0.17};
    // double qpR[] = {0.5, 0.17, 0.17, 0.17, 0.17, 0.17};
    // double qpL[] = {1, 1.54755, 0.630597, -1.06004, 1.81241, -2.4369};
    // double qpR[] = {0.5, -0.65189, -0.159904, -1.83959, -2.12555, -0.53937};
    double qpL[] = {0.999975, -0.877793, 2.07815, 1.14085, -0.926273, 0.627826};
    double qpR[] = {0.691894, -1.91336, 1.06258, -0.984366, 0.142628, 2.02722};
    double qa[] = {0, 0, 0, 0, 0, 0};

    portLogic->AddDavinciPrimitives(*T1->GetMatrix(), qpL, qa);
    portLogic->AddDavinciPrimitives(*T2->GetMatrix(), qpR, qa);
    }
}
