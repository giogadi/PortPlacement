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

#ifndef __qSlicerAutoPortPlacementFooBarWidget_h
#define __qSlicerAutoPortPlacementFooBarWidget_h

// Qt includes
#include <QWidget>

// FooBar Widgets includes
#include "qSlicerAutoPortPlacementModuleWidgetsExport.h"

class qSlicerAutoPortPlacementFooBarWidgetPrivate;

/// \ingroup Slicer_QtModules_AutoPortPlacement
class Q_SLICER_MODULE_AUTOPORTPLACEMENT_WIDGETS_EXPORT qSlicerAutoPortPlacementFooBarWidget
  : public QWidget
{
  Q_OBJECT
public:
  typedef QWidget Superclass;
  qSlicerAutoPortPlacementFooBarWidget(QWidget *parent=0);
  virtual ~qSlicerAutoPortPlacementFooBarWidget();

protected slots:

protected:
  QScopedPointer<qSlicerAutoPortPlacementFooBarWidgetPrivate> d_ptr;

private:
  Q_DECLARE_PRIVATE(qSlicerAutoPortPlacementFooBarWidget);
  Q_DISABLE_COPY(qSlicerAutoPortPlacementFooBarWidget);
};

#endif
