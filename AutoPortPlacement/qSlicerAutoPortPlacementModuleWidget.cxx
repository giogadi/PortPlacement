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

    // feasible
    // double qpL[] = {0.846626, -0.664014, -0.408133, -1.07303, -0.198602, 1.37977};
    // double qpR[] = {0.765958, -0.157138, 1.31227, 1.50681, -0.175029, 0.133613};
    // double qaL[] = {0.138726, 0.566036, -0.446022, 0.120849, -0.405235, -0.00288185};
    // double qaR[] = {-0.100714, 0.123434, 2.52511, 0.113845, 0.12373, 0.00165253};

    // feasible 2
    // double qpL[] = {0.634685, -0.125776, -0.785423, -1.34508, 0.377312, 0.754663};
    // double qpR[] = {0.609013, 0.576171, 0.522083, 2.07754, 0.549291, 0.644443};
    // double qaL[] = {0.107182, 0.677258, 0.207816, 0.191605, -0.946489, -7.88494e-11};
    // double qaR[] = {0.0342193, -0.79834, 2.37376, 0.135024, -0.59262, -6.22769e-11};

    // ut
    double qpL[] = {0.87026, 0.0410446, -2.09426, 0.589858, -0.242211, 1.61892};
    double qpR[] = {0.851171, 0.259381, 0.621712, 1.58237, 0.451428, -0.457456};
    double qaL[] = {0.0137396, 0.0202539, 0.05886, 0.114869, 0.221403, 0.011136};
    double qaR[] = {0.602573, 0.716963, 2.79583, 0.131354, 0.549648, 0.0111519};

    // initial
    // double qpL[] = {0.744029, -0.259442, -0.896604, -1.16588, 0.0253959, 0.821945};
    // double qpR[] = {0.743998, 0.258795, 0.896916, 1.16692, -0.0255285, -0.820174};
    // double qaL[] = {0.0141122, -0.024297, -0.00217849, 0.111999, 0.00570296, 0.00316633};
    // double qaR[] = {-0.0142031, -0.0244509, -3.14049, 0.111999, -0.00578722, 0.0032048};

    portLogic->AddDavinciPrimitives(*T1->GetMatrix(), qpL, qaL);
    portLogic->AddDavinciPrimitives(*T2->GetMatrix(), qpR, qaR);
    }
}
