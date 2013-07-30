/*==============================================================================

  Program: 3D Slicer

  Copyright (c) Kitware Inc.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

  This file was originally developed by Jean-Christophe Fillion-Robin, Kitware Inc.
  and was partially funded by NIH grant 3P41RR013218-12S1

==============================================================================*/

// FooBar Widgets includes
#include "qSlicerAutoPortPlacementFooBarWidget.h"
#include "ui_qSlicerAutoPortPlacementFooBarWidget.h"

//-----------------------------------------------------------------------------
/// \ingroup Slicer_QtModules_AutoPortPlacement
class qSlicerAutoPortPlacementFooBarWidgetPrivate
  : public Ui_qSlicerAutoPortPlacementFooBarWidget
{
  Q_DECLARE_PUBLIC(qSlicerAutoPortPlacementFooBarWidget);
protected:
  qSlicerAutoPortPlacementFooBarWidget* const q_ptr;

public:
  qSlicerAutoPortPlacementFooBarWidgetPrivate(
    qSlicerAutoPortPlacementFooBarWidget& object);
  virtual void setupUi(qSlicerAutoPortPlacementFooBarWidget*);
};

// --------------------------------------------------------------------------
qSlicerAutoPortPlacementFooBarWidgetPrivate
::qSlicerAutoPortPlacementFooBarWidgetPrivate(
  qSlicerAutoPortPlacementFooBarWidget& object)
  : q_ptr(&object)
{
}

// --------------------------------------------------------------------------
void qSlicerAutoPortPlacementFooBarWidgetPrivate
::setupUi(qSlicerAutoPortPlacementFooBarWidget* widget)
{
  this->Ui_qSlicerAutoPortPlacementFooBarWidget::setupUi(widget);
}

//-----------------------------------------------------------------------------
// qSlicerAutoPortPlacementFooBarWidget methods

//-----------------------------------------------------------------------------
qSlicerAutoPortPlacementFooBarWidget
::qSlicerAutoPortPlacementFooBarWidget(QWidget* parentWidget)
  : Superclass( parentWidget )
  , d_ptr( new qSlicerAutoPortPlacementFooBarWidgetPrivate(*this) )
{
  Q_D(qSlicerAutoPortPlacementFooBarWidget);
  d->setupUi(this);
}

//-----------------------------------------------------------------------------
qSlicerAutoPortPlacementFooBarWidget
::~qSlicerAutoPortPlacementFooBarWidget()
{
}
