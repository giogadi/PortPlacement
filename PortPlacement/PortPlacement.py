import os
import unittest
from __main__ import vtk, qt, ctk, slicer

#
# PortPlacement
#

class PortPlacement:
  def __init__(self, parent):
    parent.title = "Port Placement"
    parent.categories = ["IGT"]
    parent.dependencies = []
    parent.contributors = ["Andinet Enquobahrie (Kitware), Luis G. Torres (UNC)"]
    parent.helpText = string.Template("""
    The PortPlacement module assists in the planning of surgical port placement in a laparoscopic procedure. Users can specify ports using fiducial markers and the module will automatically visualize simulated surgical tools at the port locations. Users can freely pivot the simulated surgical tools about the ports, as well as automatically orient all surgical tools toward a specified surgical target. Simulated tools are represented by cylinders whose lengths and radii can be varied from tool to tool.

    See the <a href=\"$a/Documentation/Nightly/Extensions/PortPlacement\">module documentation</a> for more details.
    """).substitute({ 'a':parent.slicerWikiUrl})
    parent.acknowledgementText = """
    This work was supported by NSF GRFP Grant No. DGE-1144081 and partially supported by NIH/NIBIB Grant No. 1R43EB014074-01. Anatomical atlas volume in module and extension icons courtesy of the Surgical Planning Laboratory (<a href=\"http://www.na-mic.org/publications/item/view/1266\">Link</a>)
"""
    self.parent = parent

    # Add this test to the SelfTest module's list for discovery when the module
    # is created.  Since this module may be discovered before SelfTests itself,
    # create the list if it doesn't already exist.
    try:
      slicer.selfTests
    except AttributeError:
      slicer.selfTests = {}
    slicer.selfTests['PortPlacement'] = self.runTest

  def runTest(self):
    tester = PortPlacement()
    tester.runTest()

#
# PortPlacementWidget
#

class PortPlacementWidget:
  def __init__(self, parent = None):
    if not parent:
      self.parent = slicer.qMRMLWidget()
      self.parent.setLayout(qt.QVBoxLayout())
      self.parent.setMRMLScene(slicer.mrmlScene)
    else:
      self.parent = parent
    self.layout = self.parent.layout()
    if not parent:
      self.setup()
      self.parent.show()

  def setup(self):
    # Instantiate and connect widgets ...

    #
    # Reload and Test area
    #
    reloadCollapsibleButton = ctk.ctkCollapsibleButton()
    reloadCollapsibleButton.text = "Reload && Test"
    self.layout.addWidget(reloadCollapsibleButton)
    reloadFormLayout = qt.QFormLayout(reloadCollapsibleButton)

    # reload button
    # (use this during development, but remove it when delivering
    #  your module to users)
    self.reloadButton = qt.QPushButton("Reload")
    self.reloadButton.toolTip = "Reload this module."
    self.reloadButton.name = "PortPlacement Reload"
    reloadFormLayout.addWidget(self.reloadButton)
    self.reloadButton.connect('clicked()', self.onReload)

    # reload and test button
    # (use this during development, but remove it when delivering
    #  your module to users)
    self.reloadAndTestButton = qt.QPushButton("Reload and Test")
    self.reloadAndTestButton.toolTip = "Reload this module and then run the self tests."
    reloadFormLayout.addWidget(self.reloadAndTestButton)
    self.reloadAndTestButton.connect('clicked()', self.onReloadAndTest)

    #
    # Ports Area
    #
    portsCollapsibleButton = ctk.ctkCollapsibleButton()
    portsCollapsibleButton.text = "Surgical Ports"
    self.layout.addWidget(portsCollapsibleButton)
    portsFormLayout = qt.QFormLayout(portsCollapsibleButton)

    #
    # port fiducial list selector
    #
    self.portListSelector = slicer.qMRMLNodeComboBox()
    self.portListSelector.nodeTypes = ["vtkMRMLMarkupsFiducialNode"]
    self.portListSelector.addEnabled = False
    self.portListSelector.removeEnabled = False
    self.portListSelector.noneEnabled = True
    self.portListSelector.setMRMLScene(slicer.mrmlScene)
    self.portListSelector.setToolTip("Add surgical ports from a list of fiducials.")
    portsFormLayout.addRow("List of Surgical Ports", self.portListSelector)

    #
    # Add Port List button
    #
    self.addPortListButton = qt.QPushButton("Add Ports from List")
    self.addPortListButton.enabled = False
    portsFormLayout.addRow(self.addPortListButton)

    #
    # Port table
    #
    self.portsTable = qt.QTableView()
    self.portsTableModel = qt.QStandardItemModel()
    self.portsTable.setModel(self.portsTableModel)
    portsFormLayout.addRow(self.portsTable)

    #
    # Remove Port button
    #
    # self.removePortButton = qt.QPushButton("Remove Selected Port")
    # self.removePortButton.enabled = False
    # portsFormLayout.addRow(self.removePortButton)

    #
    # Target area
    #
    targetCollapsibleButton = ctk.ctkCollapsibleButton()
    targetCollapsibleButton.text = "Surgical Target"
    self.layout.addWidget(targetCollapsibleButton)
    targetFormLayout = qt.QFormLayout(targetCollapsibleButton)

    #
    # target selector
    #
    self.targetSelector = slicer.qMRMLNodeComboBox()
    self.targetSelector.nodeTypes = ["vtkMRMLMarkupsFiducialNode"]
    self.targetSelector.addEnabled = False
    self.targetSelector.removeEnabled = False
    self.targetSelector.noneEnabled = True
    self.targetSelector.setMRMLScene(slicer.mrmlScene)
    self.targetSelector.setToolTip("Pick the surgical target that the tools should face.")
    targetFormLayout.addRow("Surgical Target Fiducial", self.targetSelector)

    #
    # Retarget button
    #
    self.retargetButton = qt.QPushButton("Aim Tools at Target")
    self.retargetButton.toolTip = "Reset tool orientations to face target fiducial"
    self.retargetButton.enabled = False
    targetFormLayout.addRow(self.retargetButton)

    #
    # Port Tool Display Options
    #
    toolsCollapsibleButton = ctk.ctkCollapsibleButton()
    toolsCollapsibleButton.text = "Port Tool Display Options"
    self.layout.addWidget(toolsCollapsibleButton)
    toolsFormLayout = qt.QFormLayout(toolsCollapsibleButton)

    #
    # Transform sliders
    #
    # self.transformSliders = slicer.qMRMLTransformSliders()
    # self.transformSliders.TypeOfTransform = slicer.qMRMLTransformSliders.ROTATION
    # self.transformSliders.CoordinateReference = slicer.qMRMLTransformSliders.LOCAL
    # self.transformSliders.Title = 'Tool Orientation'
    # self.transformSliders.minMaxVisible = False
    # toolsFormLayout.addRow(self.transformSliders)

    #
    # radius spin box
    #
    self.radiusSpinBox = qt.QDoubleSpinBox()
    self.radiusSpinBox.setMinimum(0.0)
    self.radiusSpinBox.setMaximum(10.0)
    self.radiusSpinBox.setValue(2.0)
    toolsFormLayout.addRow("Tool radius", self.radiusSpinBox)

    #
    # length spin box
    #
    self.lengthSpinBox = qt.QDoubleSpinBox()
    self.lengthSpinBox.setMinimum(0.0)
    self.lengthSpinBox.setMaximum(250.0)
    self.lengthSpinBox.setValue(150.0)
    toolsFormLayout.addRow("Tool length", self.lengthSpinBox)

    # connections
    self.portListSelector.connect('currentNodeChanged(bool)', self.onPortListSelectorChanged)
    self.targetSelector.connect('currentNodeChanged(bool)', self.onTargetSelectorChanged)
    self.addPortListButton.connect('clicked(bool)', self.onAddPortListButton)
    # self.removePortButton.connect('clicked(bool)', self.onRemovePortButton)
    self.retargetButton.connect('clicked(bool)', self.onRetargetButton)
    # self.radiusSpinBox.connect('valueChanged(double)', self.onToolShapeChanged)
    # self.lengthSpinBox.connect('valueChanged(double)', self.onToolShapeChanged)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Add observer to scene for removal of fiducial nodes
    self.sceneObserverTag = slicer.mrmlScene.AddObserver(slicer.mrmlScene.NodeRemovedEvent,
                                                         self.onNodeRemoved)

    # instantiate port placement module logic
    self.logic = PortPlacementLogic()

  def __del__(self):
    slicer.mrmlScene.RemoveObserver(self.sceneObserverTag)

  def onPortListSelectorChanged(self, valid):
    self.addPortListButton.enabled = valid

  def onTargetSelectorChanged(self, valid):
    self.retargetButton.enabled = valid

  def onAddPortListButton(self):
    self.logic.setMarkupsNode(self.portListSelector.currentNode(),
                              self.radiusSpinBox.value,
                              self.lengthSpinBox.value)
    # (portFiducials, portVisibility) = self.logic.addTools(self.portListSelector.currentNode(),
    #                                                       self.radiusSpinBox.value,
    #                                                       self.lengthSpinBox.value)
    # self.updateTable(portFiducials, portVisibility)

  # def onRemovePortButton(self):
  #   # Problem: currentIndex() can be valid even when there is no
  #   # visible selection. This will have to do for now while I figure
  #   # out how to making Python bindings to QList<QModelIndex> (ugh)
  #   #
  #   # TODO: instead of using index in table, use our new item ->
  #   # fiducial mapping
  #   index = self.portsTable.selectionModel().currentIndex()
  #   (portFiducials, portVisibility) = self.logic.removeToolByIndex(index.row())
  #   self.updateTable(portFiducials, portVisibility)

  # def onTableItemChanged(self,item):
  #   # Get fiducial corresponding to the table item
  #   fiducial = self.itemFiducialMap[item]

  #   # If the item is checkable, it's a visibility column. If not, it's
  #   # the name column
  #   if item.isCheckable():
  #     if item.checkState() == 2: # checked enum
  #       self.logic.makeToolVisible(fiducial)
  #     else:
  #       self.logic.makeToolInvisible(fiducial)

  #   else: # name column
  #     fiducial.SetName(item.text())

  def onNodeRemoved(self, scene, event):
    self.logic.onNodeRemoved()
    # (portFiducials, portVisibility) = self.logic.onNodeRemoved()
    # self.updateTable(portFiducials, portVisibility)

  def onRetargetButton(self):
    self.logic.retargetTools(self.targetSelector.currentNode())

  # def onToolShapeChanged(self):
  #     toolIdx = self.portsTable.selectionModel().currentIndex().row()
  #     if toolIdx < self.portsTableModel.rowCount:
  #       self.logic.setToolShapeByIndex(toolIdx, self.radiusSpinBox.value, self.lengthSpinBox.value)

  # def onCurrentToolChanged(self, newIndex, prevIndex):
  #   self.removePortButton.enabled = True
  #   self.transformSliders.setMRMLTransformNode(self.logic.getToolTransformByIndex(newIndex.row()))
  #   self.radiusSpinBox.setValue(self.logic.getToolRadiusByIndex(newIndex.row()))
  #   self.lengthSpinBox.setValue(self.logic.getToolLengthByIndex(newIndex.row()))

  # def updateTable(self, portFiducials, portVisibility):
  #   self.itemFiducialMap = {}
  #   self.portsTableModel = qt.QStandardItemModel()
  #   self.portsTable.setModel(self.portsTableModel)
  #   for (row,(fid, vis)) in enumerate(zip(portFiducials,portVisibility)):
  #     item = qt.QStandardItem()
  #     item.setText(fid.GetName())
  #     self.portsTableModel.setItem(row,0,item)
  #     self.itemFiducialMap[item] = fid

  #     item = qt.QStandardItem()
  #     item.setText('Visible')
  #     item.setCheckable(True)
  #     item.setCheckState(vis)
  #     self.portsTableModel.setItem(row,1,item)
  #     self.itemFiducialMap[item] = fid
  #   self.portsTableModel.setHeaderData(0,1,"Port Fiducial Name")
  #   self.portsTableModel.setHeaderData(1,1," ")
  #   self.portsTable.setColumnWidth(0,15*len("Port Fiducial Name"))

  #   self.portsTableModel.itemChanged.connect(self.onTableItemChanged)
  #   self.transformSliders.reset()
  #   self.removePortButton.enabled = False
  #   self.transformSliders.setMRMLTransformNode(None)
  #   self.portsTable.selectionModel().currentRowChanged.connect(self.onCurrentToolChanged)

  def onReload(self,moduleName="PortPlacement"):
    """Generic reload method for any scripted module.
    ModuleWizard will subsitute correct default moduleName.
    """
    import imp, sys, os, slicer

    widgetName = moduleName + "Widget"

    # reload the source code
    # - set source file path
    # - load the module to the global space
    filePath = eval('slicer.modules.%s.path' % moduleName.lower())
    p = os.path.dirname(filePath)
    if not sys.path.__contains__(p):
      sys.path.insert(0,p)
    fp = open(filePath, "r")
    globals()[moduleName] = imp.load_module(
        moduleName, fp, filePath, ('.py', 'r', imp.PY_SOURCE))
    fp.close()

    # rebuild the widget
    # - find and hide the existing widget
    # - create a new widget in the existing parent
    parent = slicer.util.findChildren(name='%s Reload' % moduleName)[0].parent().parent()
    for child in parent.children():
      try:
        child.hide()
      except AttributeError:
        pass
    # Remove spacer items
    item = parent.layout().itemAt(0)
    while item:
      parent.layout().removeItem(item)
      item = parent.layout().itemAt(0)
    # create new widget inside existing parent
    globals()[widgetName.lower()] = eval(
        'globals()["%s"].%s(parent)' % (moduleName, widgetName))
    globals()[widgetName.lower()].setup()

  def onReloadAndTest(self,moduleName="PortPlacement"):
    try:
      self.onReload()
      evalString = 'globals()["%s"].%sTest()' % (moduleName, moduleName)
      tester = eval(evalString)
      tester.runTest()
    except Exception, e:
      import traceback
      traceback.print_exc()
      qt.QMessageBox.warning(slicer.util.mainWindow(),
          "Reload and Test", 'Exception!\n\n' + str(e) + "\n\nSee Python Console for Stack Trace")


#
# PortPlacementLogic
#

class PortPlacementLogic:
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget
  """

  class Tool:
    def __init__(self, markupID, radius, length, position):
      self.markupID = markupID

      # Create a vtkCylinderSource for rendering the tool
      self.toolSrc = vtk.vtkCylinderSource()
      self.toolSrc.SetRadius(radius)
      self.toolSrc.SetHeight(length)

      # Create tool model using toolSrc
      self.modelNode = slicer.vtkMRMLModelNode()
      self.modelNode.SetName(slicer.mrmlScene.GenerateUniqueName("Tool"))
      polyData = self.toolSrc.GetOutput()
      self.modelNode.SetAndObservePolyData(polyData)
      slicer.mrmlScene.AddNode(self.modelNode)

      # Create a DisplayNode for this node (so that it can control its
      # own visibility) and have modelNode observe it
      # and then our model display node
      self.modelDisplay = slicer.vtkMRMLModelDisplayNode()
      self.modelDisplay.SetColor(0,1,1) # cyan
      slicer.mrmlScene.AddNode(self.modelDisplay)
      self.modelNode.SetAndObserveDisplayNodeID(self.modelDisplay.GetID())

      # Create a transform for our tool model that will scale it to
      # the proper radius, length, and position. Leave the orientation
      # at the default.
      self.transformNode = slicer.vtkMRMLLinearTransformNode()
      self.transformNode.SetName(slicer.mrmlScene.GenerateUniqueName("Transform"))
      self.__updateTransformNode__(position)
      slicer.mrmlScene.AddNode(self.transformNode)

      # apply the new transform to our tool
      self.modelNode.SetAndObserveTransformNodeID(self.transformNode.GetID())

      # Set tool to be visible
      self.visible = True

    def updateShape(self, radius, length):
      self.toolSrc.SetRadius(radius)
      self.toolSrc.SetHeight(length)
      self.toolSrc.Update()

    def updatePosition(self, position):
      self.__updateTransformNode__(position)

    def makeVisible(self):
      self.modelDisplay.VisibilityOn()
      self.visible = True

    def makeInvisible(self):
      self.modelDisplay.VisibilityOff()
      self.visible = False

    def __del__(self):
      slicer.mrmlScene.RemoveNode(self.modelNode)
      slicer.mrmlScene.RemoveNode(self.modelDisplay)
      slicer.mrmlScene.RemoveNode(self.transformNode)

    def __updateTransformNode__(self, p):
      # First get the transformNode's current transform it, invert it
      # to get back to identity, then set the new transform
      #
      # TODO: figure out if there's a better way to do this, and if
      # not, report an issue about this.
      mat = vtk.vtkMatrix4x4()
      self.transformNode.GetMatrixTransformToWorld(mat)

      t = vtk.vtkTransform()
      t.SetMatrix(mat)
      orientWXYZ = [0,0,0,0]
      t.GetOrientationWXYZ(orientWXYZ)

      t.Inverse()
      self.transformNode.ApplyTransformMatrix(t.GetMatrix())

      t.Identity()
      t.Translate(p)
      t.RotateWXYZ(orientWXYZ[0], orientWXYZ[1:])
      self.transformNode.ApplyTransformMatrix(t.GetMatrix())

  def __init__(self):
    self.markupsNode = None
    self.toolList = []

  def __del__(self):
    if not (self.markupsNode is None):
      self.MarkupsNode.RemoveObserver(self.nodeObserverTag)

  def clearMarkupsNode(self):
    # Clear the current tool list
    self.toolList = []

    # If we were observing another node, remove first remove ourself
    # as an observer
    if not (self.markupsNode is None):
      self.markupsNode.RemoveObserver(self.pointModifiedObserverTag)
      self.markupsNode.RemoveObserver(self.markupRemovedObserverTag)

  def setMarkupsNode(self, node, radius, length):
    self.clearMarkupsNode()

    self.markupsNode = node
    if not (node is None):
      self.pointModifiedObserverTag = node.AddObserver(slicer.vtkMRMLMarkupsNode.PointModifiedEvent, self.onMarkupModified)
      self.markupRemovedObserverTag = node.AddObserver(slicer.vtkMRMLMarkupsNode.MarkupRemovedEvent, self.onMarkupRemoved)

      # Add tools associated with the markups in markupsNode (if type is
      # vtkMRMLMarkupsFiducial)
      for i in range(node.GetNumberOfFiducials()):
        p = [0,0,0]
        node.GetNthFiducialPosition(i, p)
        self.toolList.append(self.Tool(node.GetNthMarkupID(i), radius, length, p))

  def setToolShapeByIndex(self, idx, radius, length):
    self.toolList[idx].updateShape(radius, length)

  def makeToolVisible(self, idx):
    self.toolList[idx].makeVisible()

  def makeToolInvisible(self, idx):
    self.toolList[idx].makeInvisible()

  def getToolTransformByIndex(self, idx):
    return self.toolList[idx].transformNode

  def getToolRadiusByIndex(self, idx):
    return self.toolList[idx].radius

  def getToolLengthByIndex(self, idx):
    return self.toolList[idx].length

  # If the node we're observing has been modified, just iterate
  # through and update all the tools since for now we don't know how
  # to distinguish which markup got modified.
  def onMarkupModified(self, node, event):
    for i in range(node.GetNumberOfFiducials()):
      p = [0,0,0]
      node.GetNthFiducialPosition(i, p)
      self.toolList[i].updatePosition(p)

  def onNodeRemoved(self):
    if (self.markupsNode is None):
      return

    if not slicer.mrmlScene.IsNodePresent(self.markupsNode):
      self.clearMarkupsNode()

  def onMarkupRemoved(self, node, event):
    if self.markupsNode is None:
      return

    markupsSet = set()
    for i in range(self.markupsNode.GetNumberOfFiducials()):
      markupsSet.add(self.markupsNode.GetNthMarkupID(i))

    prevLen = len(self.toolList)
    for i in range(len(self.toolList)):
      if not (self.toolList[i].markupID in markupsSet):
        del self.toolList[i]
        break
    if prevLen == len(self.toolList):
      print("One of our nodes wasn't removed...")

  # These flags are a horrible, HORRIBLE result of not being able to
  # inspect what kind of node the scene just removed (because it
  # hasn't yet been wrapped in Python). If the flags weren't there,
  # what would happen is that upon removal of the tool model, this
  # function would get called AGAIN, and so-on for infinite
  # recursion. So, with the flags, we guarantee that the body of this
  # function only gets called for the fiducial node removal.
  #
  # We also want to look for removal of the currently active node
  # hierarchy.
  # def onNodeRemoved(self):
  #   if self.flag:
  #     self.flag = False
  #     # This seems like a really wasteful solution, but for now it's
  #     # necessary because apparently the calldata representing the
  #     # removed fiducial has not yet been wrapped in Python.
  #     for fiducial in self.fiducialToolMap.keys():
  #       if not slicer.mrmlScene.IsNodePresent(fiducial):
  #         self.removeTool(fiducial)
  #     self.flag = True
  #   return self.getPortsAndVisibility()

  # def addTool(self, fiducialNode, radius, length):
  #   # Make sure this annotation is indeed a fiducial marker, and don't
  #   # add the same fiducial twice
  #   if fiducialNode.GetClassName() == 'vtkMRMLAnnotationFiducialNode' and \
  #         not fiducialNode in self.fiducialToolMap:
  #     # create the tool model
  #     toolModel = slicer.vtkMRMLModelNode()
  #     toolModel.SetName(slicer.mrmlScene.GenerateUniqueName("Tool"))
  #     toolSrc = vtk.vtkCylinderSource()
  #     toolSrc.SetRadius(radius)
  #     toolSrc.SetHeight(length)
  #     polyData = toolSrc.GetOutput()
  #     toolModel.SetAndObservePolyData(polyData)
  #     slicer.mrmlScene.AddNode(toolModel)

  #     # and then our model display node
  #     modelDisplay = slicer.vtkMRMLModelDisplayNode()
  #     modelDisplay.SetColor(0,1,1) # cyan
  #     slicer.mrmlScene.AddNode(modelDisplay)
  #     toolModel.SetAndObserveDisplayNodeID(modelDisplay.GetID())

  #     # and our transform node to correspond with the fiducial's position
  #     fidPos = [0,0,0]
  #     fiducialNode.GetFiducialCoordinates(fidPos)
  #     t = vtk.vtkTransform()
  #     t.Translate(fidPos)
  #     transformNode = slicer.vtkMRMLLinearTransformNode()
  #     transformNode.ApplyTransformMatrix(t.GetMatrix())
  #     transformNode.SetName(slicer.mrmlScene.GenerateUniqueName("Transform"))
  #     slicer.mrmlScene.AddNode(transformNode)

  #     # apply the new transform to our tool
  #     toolModel.SetAndObserveTransformNodeID(transformNode.GetID())

  #     # add an observer to this fiducial node to monitor changes in its position
  #     tag = fiducialNode.AddObserver('ModifiedEvent', self.onFiducialModified)

  #     # add the model, model display, and transform to our fiducialToolMap
  #     self.fiducialToolMap[fiducialNode] = self.Tool(toolSrc,
  #                                                    toolModel,
  #                                                    modelDisplay,
  #                                                    True,
  #                                                    transformNode,
  #                                                    tag)
  #   return self.getPortsAndVisibility()

  # def addTools(self, annotationHierarchy, radius, length):
  #   fiducialList = self.fromHierarchyToFiducialList(annotationHierarchy)
  #   for f in fiducialList:
  #     self.addTool(f, radius, length)
  #   return self.getPortsAndVisibility()

  # def getPortsAndVisibility(self):
  #   return (self.fiducialToolMap.keys(),
  #           [self.fiducialToolMap[f].visible for f in self.fiducialToolMap])

  # def fromHierarchyToFiducialList(self, annotationHierarchy):
  #   collection = vtk.vtkCollection()
  #   annotationHierarchy.GetChildrenDisplayableNodes(collection)
  #   portFiducialList = [collection.GetItemAsObject(i) for i in
  #                       xrange(collection.GetNumberOfItems())]
  #   return [fid for fid in portFiducialList
  #           if fid.GetClassName() == "vtkMRMLAnnotationFiducialNode"]

  def retargetTools(self, targetFid):

    if self.markupsNode is None or \
       targetFid.GetNumberOfFiducials() < 1 or \
       targetFid == self.markupsNode:
      return

    import numpy
    import numpy.linalg
    import math
    # Get target coordinates
    targetLoc = [0,0,0]
    targetFid.GetNthFiducialPosition(0, targetLoc)
    targetLoc = numpy.array(targetLoc)

    # Iterate through port tool map and retarget their associated tools
    for i in range(self.markupsNode.GetNumberOfFiducials()):
      portLocList = [0,0,0]
      self.markupsNode.GetNthFiducialPosition(i, portLocList)
      portLoc = numpy.array(portLocList)

      # Assume tools get drawn aligned with the global y-axis. We
      # want to transform the tool to be oriented along the port
      # toward the target point. We begin by finding the axis to
      # rotate the tool by, which is the cross product of the
      # global y and the target vector
      targetVec = targetLoc - portLoc

      if numpy.linalg.norm(targetVec) <= 0.0000001:
        continue

      targetVec = targetVec / numpy.linalg.norm(targetVec)
      rotAxis = numpy.cross(numpy.array([0,1,0]), targetVec)
      rotAxis = rotAxis / numpy.linalg.norm(rotAxis)

      # get rotation angle and normalize rotation axis
      angle = math.acos(numpy.dot(numpy.array([0,1,0]), targetVec)) * 180. / math.pi

      # generate our transform
      #
      # We do this by first undo-ing the tool's current transform and
      # then applying the new one.
      mat = vtk.vtkMatrix4x4()
      self.toolList[i].transformNode.GetMatrixTransformToWorld(mat)
      t = vtk.vtkTransform()
      t.SetMatrix(mat)

      t.Inverse()
      self.toolList[i].transformNode.ApplyTransformMatrix(t.GetMatrix())

      t.Identity()
      t.Translate(portLocList)
      t.RotateWXYZ(angle, rotAxis.tolist())
      self.toolList[i].transformNode.ApplyTransformMatrix(t.GetMatrix())


class PortPlacementTest(unittest.TestCase):
  """
  This is the test case for your scripted module.
  """

  def delayDisplay(self,message,msec=1000):
    """This utility method displays a small dialog and waits.
    This does two things: 1) it lets the event loop catch up
    to the state of the test so that rendering and widget updates
    have all taken place before the test continues and 2) it
    shows the user/developer/tester the state of the test
    so that we'll know when it breaks.
    """
    print(message)
    self.info = qt.QDialog()
    self.infoLayout = qt.QVBoxLayout()
    self.info.setLayout(self.infoLayout)
    self.label = qt.QLabel(message,self.info)
    self.infoLayout.addWidget(self.label)
    qt.QTimer.singleShot(msec, self.info.close)
    self.info.exec_()

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    # self.test_PortPlacement1()

  # def test_PortPlacement1(self):
  #   import numpy
  #   import numpy.linalg
  #   import random

  #   self.delayDisplay("Starting test...")

  #   m = slicer.util.mainWindow()
  #   m.moduleSelector().selectModule('PortPlacement')
  #   self.widget = slicer.modules.PortPlacementWidget

  #   logic = self.widget.logic

  #   annotationsLogic = slicer.util.getModule('Annotations').logic()
  #   f = slicer.vtkMRMLAnnotationFiducialNode()
  #   f.SetFiducialCoordinates(*[random.uniform(-100.,100.) for j in range(3)])
  #   f.Initialize(slicer.mrmlScene)

  #   self.assertTrue(not self.widget.addPortButton.enabled)
  #   self.widget.portFiducialSelector.setCurrentNode(f)
  #   self.assertTrue(self.widget.addPortButton.enabled)
  #   self.widget.onAddPortButton()
  #   self.widget.onAddPortButton()

  #   # check table
  #   self.assertTrue(self.widget.portsTableModel.rowCount() == 1)
  #   self.assertTrue(self.widget.portsTableModel.item(0).text() == f.GetName())

  #   # check tool transform
  #   center = [0,0,0,1]
  #   toolPosition = [0,0,0,1]
  #   logic.fiducialToolMap[f].model.TransformPointToWorld(center, toolPosition)
  #   fidCoords = [0,0,0]
  #   f.GetFiducialCoordinates(fidCoords)
  #   diff = numpy.array(toolPosition)[0:3] - numpy.array(fidCoords)
  #   self.assertTrue(numpy.dot(diff,diff) < 1e-10)

  #   index = self.widget.portsTableModel.index(0,0)
  #   self.assertTrue(not self.widget.removePortButton.enabled)
  #   self.widget.portsTable.selectionModel().setCurrentIndex(index, qt.QItemSelectionModel.Current)
  #   self.assertTrue(self.widget.removePortButton.enabled)
  #   self.widget.onRemovePortButton()
  #   self.assertTrue(self.widget.portsTableModel.rowCount() == 0)
  #   self.assertTrue(len(logic.fiducialToolMap) == 0)

  #   # add a list of ports
  #   annotationsLogic.SetActiveHierarchyNodeID(annotationsLogic.GetTopLevelHierarchyNodeID())
  #   annotationsLogic.AddHierarchy()
  #   portsHierarchy = annotationsLogic.GetActiveHierarchyNode()

  #   # Add some ports to our ports hierarchy
  #   portFiducialNodes = []
  #   numPorts = 4
  #   for i in range(numPorts):
  #     fiducialNode = slicer.vtkMRMLAnnotationFiducialNode()
  #     fiducialNode.SetFiducialCoordinates(*[random.uniform(-100.,100.) for j in range(3)])
  #     fiducialNode.Initialize(slicer.mrmlScene)
  #     portFiducialNodes.append(fiducialNode)

  #   self.assertTrue(not self.widget.addPortListButton.enabled)
  #   self.widget.portListSelector.setCurrentNode(portsHierarchy)
  #   self.assertTrue(self.widget.portListSelector.enabled)
  #   self.widget.onAddPortListButton()
  #   self.widget.onAddPortListButton()

  #   # check new fiducials
  #   self.assertTrue(self.widget.portsTableModel.rowCount() == numPorts)
  #   keys = logic.fiducialToolMap.keys()
  #   self.assertTrue(len(keys) == numPorts)
  #   for (i,key) in enumerate(keys):
  #     self.assertTrue(key in portFiducialNodes)
  #     self.assertTrue(self.widget.portsTableModel.item(i).text() == key.GetName())

  #     # check tool transform
  #     logic.fiducialToolMap[key].model.TransformPointToWorld(center, toolPosition)
  #     key.GetFiducialCoordinates(fidCoords)
  #     diff = numpy.array(toolPosition)[0:3] - numpy.array(fidCoords)
  #     self.assertTrue(numpy.dot(diff,diff) < 1e-10)

  #   # check retargeting
  #   targetFiducial = slicer.vtkMRMLAnnotationFiducialNode()
  #   targetFiducial.SetFiducialCoordinates(*[random.uniform(-100.,100.) for i in range(3)])
  #   annotationsLogic.SetActiveHierarchyNodeID(annotationsLogic.GetTopLevelHierarchyNodeID())
  #   targetFiducial.Initialize(slicer.mrmlScene)
  #   annotationsLogic.SetActiveHierarchyNodeID(portsHierarchy.GetID())

  #   self.assertTrue(not self.widget.retargetButton.enabled)
  #   self.widget.targetSelector.setCurrentNode(targetFiducial)
  #   self.assertTrue(self.widget.retargetButton.enabled)
  #   self.widget.onRetargetButton()

  #   # check retargeting by verifying that tools' positions are
  #   # unchanged and that their y-axes are oriented toward point
  #   for key in keys:
  #     logic.fiducialToolMap[key].model.TransformPointToWorld(center, toolPosition)
  #     key.GetFiducialCoordinates(fidCoords)
  #     diff = numpy.array(toolPosition)[0:3] - numpy.array(fidCoords)
  #     self.assertTrue(numpy.dot(diff,diff) < 1e-10)

  #     targetWorld = [0,0,0]
  #     targetFiducial.GetFiducialCoordinates(targetWorld)
  #     targetWorld = targetWorld + [1]
  #     targetLocal = [0,0,0,1]
  #     logic.fiducialToolMap[key].model.TransformPointFromWorld(targetWorld, targetLocal)

  #     targetLocal = numpy.array(targetLocal)[0:3]
  #     targetLocal = targetLocal / numpy.linalg.norm(targetLocal)

  #     # target local should be the unit y-axis (e2)
  #     yAxis = numpy.array([0.,1.,0.])
  #     diff = yAxis - targetLocal
  #     self.assertTrue(numpy.dot(diff,diff) < 1e-10)

  #   # dunno how to test sliders...

  #   self.delayDisplay("Test passed!")
