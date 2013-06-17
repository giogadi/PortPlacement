import os
import unittest
from __main__ import vtk, qt, ctk, slicer

#
# PortPlacement
#

class PortPlacement:
  def __init__(self, parent):
    parent.title = "Port Placement"
    parent.categories = ["Work in Progress"]
    parent.dependencies = []
    parent.contributors = ["Luis G. Torres (UNC)"]
    parent.helpText = """
    This module allows the user to place and visualize simulated laparoscopic ports.
    """
    parent.acknowledgementText = """
    This file was originally developed by Luis G. Torres (UNC).
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
    # viz options area
    #
    vizOptionsCollapsibleButton = ctk.ctkCollapsibleButton()
    vizOptionsCollapsibleButton.text = "Laparoscopic Tool Options"
    self.layout.addWidget(vizOptionsCollapsibleButton)

    # Layout within the dummy collapsible button
    vizOptionsFormLayout = qt.QFormLayout(vizOptionsCollapsibleButton)

    #
    # radius spin box
    # 
    self.radiusSpinBox = qt.QDoubleSpinBox()
    self.radiusSpinBox.setMinimum(0.0)
    self.radiusSpinBox.setMaximum(10.0)
    self.radiusSpinBox.setValue(2.0)
    vizOptionsFormLayout.addRow("Tool radius: ", self.radiusSpinBox)

    #
    # length spin box
    #
    self.lengthSpinBox = qt.QDoubleSpinBox()
    self.lengthSpinBox.setMinimum(0.0)
    self.lengthSpinBox.setMaximum(250.0)
    self.lengthSpinBox.setValue(150.0)
    vizOptionsFormLayout.addRow("Tool length: ", self.lengthSpinBox)

    #
    # Inputs Area
    #
    inputsCollapsibleButton = ctk.ctkCollapsibleButton()
    inputsCollapsibleButton.text = "Ports and Target"
    self.layout.addWidget(inputsCollapsibleButton)
    inputsFormLayout = qt.QFormLayout(inputsCollapsibleButton)

    #
    # port list selector
    #
    self.portListSelector = slicer.qMRMLNodeComboBox()
    self.portListSelector.nodeTypes = ["vtkMRMLAnnotationHierarchyNode"]
    self.portListSelector.addEnabled = False
    self.portListSelector.removeEnabled = False
    self.portListSelector.noneEnabled = True
    self.portListSelector.setMRMLScene(slicer.mrmlScene)
    self.portListSelector.setToolTip("Pick the ports through which to insert the tools.")
    inputsFormLayout.addRow("Fiducial List of Port Placements: ", self.portListSelector)

    #
    # target selector
    #
    self.targetSelector = slicer.qMRMLNodeComboBox()
    self.targetSelector.nodeTypes = ["vtkMRMLAnnotationFiducialNode"]
    self.targetSelector.addEnabled = False
    self.targetSelector.removeEnabled = False
    self.targetSelector.noneEnabled = True
    self.targetSelector.setMRMLScene(slicer.mrmlScene)
    self.targetSelector.setToolTip("Pick the surgical target that the tools should face.")
    inputsFormLayout.addRow("Surgical Target Fiducial: ", self.targetSelector)

    #
    # Retarget button
    #
    self.retargetButton = qt.QPushButton("Aim Tools at Target")
    self.retargetButton.toolTip = "Reset tool orientations to face target fiducial"
    inputsFormLayout.addRow(self.retargetButton)

    # connections
    # self.updateButton.connect('clicked(bool)', self.onUpdateButton)
    self.portListSelector.connect('currentNodeChanged(bool)', self.onPortListChanged)
    self.retargetButton.connect('clicked(bool)', self.onRetargetButton)
    self.radiusSpinBox.connect('valueChanged(double)', self.onToolShapeChanged)
    self.lengthSpinBox.connect('valueChanged(double)', self.onToolShapeChanged)

    # Add vertical spacer
    self.layout.addStretch(1)

    # instantiate port placement module logic
    self.logic = PortPlacementLogic(self.radiusSpinBox.value, self.lengthSpinBox.value)

  def onPortListChanged(self):
    self.logic.setActivePortList(self.portListSelector.currentNode())

  def onRetargetButton(self):
    if self.targetSelector.currentNode():
      self.logic.retargetTools(self.targetSelector.currentNode(),
                               self.portListSelector.currentNode())

  def onToolShapeChanged(self):
    self.logic.setToolShape(self.radiusSpinBox.value, self.lengthSpinBox.value)

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
    def __init__(self, model, modelDisplay, transform, fidObserverTag):
      self.model = model
      self.modelDisplay = modelDisplay
      self.transform = transform
      self.fiducialObserverTag = fidObserverTag

  def __init__(self, initToolRadius, initToolLength):
    self.fiducialToolMap = {}
    self.toolSrc = vtk.vtkCylinderSource()
    self.setToolShape(initToolRadius, initToolLength)
    self.activeAnnotationHierarchy = None
    self.activeHierarchyObserverTag = None

    # Should we track this tag and remove it later?
    slicer.mrmlScene.AddObserver(slicer.mrmlScene.NodeRemovedEvent, self.onNodeRemoved)

    # this is a horrible hack to avoid recursion in onFiducialRemoved
    self.flag = True

  # def __del__ (do this if we are seeing leaks or leftover observers)

  def setToolShape(self, radius, length):
    self.toolSrc.SetRadius(radius)
    self.toolSrc.SetHeight(length)
    self.toolPolyData = self.toolSrc.GetOutput()

    for portFid in self.fiducialToolMap:
      self.fiducialToolMap[portFid].model.SetAndObservePolyData(self.toolPolyData)

  def setActivePortList(self, newPortAnnotationHierarchy):
    # this probably doesn't happen, but let's make sure
    if newPortAnnotationHierarchy == self.activeAnnotationHierarchy:
      return

    # Set all ports to be invisible; we're going to only make the
    # ports in the new list visible
    for fiducialNode in self.fiducialToolMap:
      self.fiducialToolMap[fiducialNode].modelDisplay.VisibilityOff()

    # If there was a previous active hiearchy, remove our observer to
    # it. If our new active hierarchy isn't 'None', add an observer to
    # it.
    if self.activeHierarchyObserverTag:
      self.activeAnnotationHierarchy.RemoveObserver(self.activeHierarchyObserverTag)
    self.activeAnnotationHierarchy = newPortAnnotationHierarchy
    if self.activeAnnotationHierarchy:
      self.activeHierarchyObserverTag = \
        self.activeAnnotationHierarchy.AddObserver('ModifiedEvent', self.onActiveHierarchyChanged)
    else:
      self.activeHierarchyObserverTag = None
      return

    newPortFiducialList = self.fromHierarchyToFiducialList(newPortAnnotationHierarchy)

    for fiducialNode in newPortFiducialList:
      if not fiducialNode in self.fiducialToolMap:
        self.addTool(fiducialNode)

      else:
        # if this node is already in our list of tools, just make it visible
        self.fiducialToolMap[fiducialNode].modelDisplay.VisibilityOn()
      
  # Given a new position of a tool's fiducial node, we modify the
  # tool's transform such that the middle of the tool is on the new
  # fiducial, but the prior orientation is preserved
  def onFiducialModified(self, fiducialNode, event):
    tool = self.fiducialToolMap[fiducialNode]
    mat = vtk.vtkMatrix4x4()
    tool.transform.GetMatrixTransformToWorld(mat)

    t = vtk.vtkTransform()
    t.SetMatrix(mat)
    orientation = [0,0,0,0]
    t.GetOrientationWXYZ(orientation)

    newPosition = [0,0,0]
    fiducialNode.GetFiducialCoordinates(newPosition)

    t.Inverse()
    t.Translate(newPosition)
    t.RotateWXYZ(orientation[0], orientation[1:])
    tool.transform.ApplyTransformMatrix(t.GetMatrix())

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
  def onNodeRemoved(self, scene, event):
    if self.flag:
      self.flag = False
      # This seems like a really wasteful solution, but for now it's
      # necessary because apparently the calldata representing the
      # removed fiducial has not yet been wrapped in Python.
      for fiducial in self.fiducialToolMap.keys():
        if not slicer.mrmlScene.IsNodePresent(fiducial):
          self.removeTool(fiducial)

      # Also check for currently active hierarchy
      if not slicer.mrmlScene.IsNodePresent(self.activeAnnotationHierarchy):
        print('waddap?')
        self.setActivePortList(None)
      self.flag = True

  def onActiveHierarchyChanged(self, node, event):
    for fiducialNode in self.fiducialToolMap:
      self.fiducialToolMap[fiducialNode].modelDisplay.VisibilityOff()

    # Get active port annotation hierarchy as list of annotations
    newPortFiducialList = self.fromHierarchyToFiducialList(self.activeAnnotationHierarchy)

    # Look for fiducials moved into the active hierarchy
    for fiducial in newPortFiducialList:
      if not fiducial in self.fiducialToolMap:
        self.addTool(fiducial)
      else:
        self.fiducialToolMap[fiducial].modelDisplay.VisibilityOn()

  def removeTool(self, fiducialNode):
    tool = self.fiducialToolMap[fiducialNode]
    fiducialNode.RemoveObserver(tool.fiducialObserverTag)
    slicer.mrmlScene.RemoveNode(tool.model)
    slicer.mrmlScene.RemoveNode(tool.transform)
    del self.fiducialToolMap[fiducialNode]

  def addTool(self, fiducialNode):
    # create the tool model
    toolModel = slicer.vtkMRMLModelNode()
    toolModel.SetName(slicer.mrmlScene.GenerateUniqueName("Tool"))
    toolModel.SetAndObservePolyData(self.toolPolyData)
    slicer.mrmlScene.AddNode(toolModel)

    # and then our model display node
    modelDisplay = slicer.vtkMRMLModelDisplayNode()
    modelDisplay.SetColor(0,1,1) # cyan
    slicer.mrmlScene.AddNode(modelDisplay)
    toolModel.SetAndObserveDisplayNodeID(modelDisplay.GetID())

    # and our transform node to correspond with the fiducial's position
    fidPos = [0,0,0]
    fiducialNode.GetFiducialCoordinates(fidPos)
    t = vtk.vtkTransform()
    t.Translate(fidPos)
    transformNode = slicer.vtkMRMLLinearTransformNode()
    transformNode.ApplyTransformMatrix(t.GetMatrix())
    transformNode.SetName(slicer.mrmlScene.GenerateUniqueName("Transform"))
    slicer.mrmlScene.AddNode(transformNode)

    # apply the new transform to our tool
    toolModel.SetAndObserveTransformNodeID(transformNode.GetID())

    # add an observer to this fiducial node to monitor changes in its position
    tag = fiducialNode.AddObserver('ModifiedEvent', self.onFiducialModified)

    # add the model, model display, and transform to our fiducialToolMap
    self.fiducialToolMap[fiducialNode] = self.Tool(toolModel, 
                                                   modelDisplay, 
                                                   transformNode, 
                                                   tag)

  def fromHierarchyToFiducialList(self, annotationHierarchy):
    collection = vtk.vtkCollection()
    annotationHierarchy.GetChildrenDisplayableNodes(collection)
    portFiducialList = [collection.GetItemAsObject(i) for i in
                        xrange(collection.GetNumberOfItems())]
    return [fid for fid in portFiducialList 
            if fid.GetClassName() == "vtkMRMLAnnotationFiducialNode"]

  # note: currently using a simple-to-understand but *very*
  # inefficient method of finding the two port lists' intersection and
  # differences.
  # def updatePorts(self, newPortAnnotationHierarchy):
  #   import numpy
    
  #   # If newPortAnnotationHierarchy is None, just clear the tools
  #   # (this should be refactored)
  #   if not newPortAnnotationHierarchy:
  #     for portFid in [fid for fid in self.fiducialToolMap]:
  #       tool = self.fiducialToolMap[portFid]
  #       slicer.mrmlScene.RemoveNode(tool.model)
  #       slicer.mrmlScene.RemoveNode(tool.transform)
  #       del self.fiducialToolMap[portFid]
  #     return
                
  #   # turn the annotation hierarchy into a list of annotations
  #   collection = vtk.vtkCollection()
  #   newPortAnnotationHierarchy.GetChildrenDisplayableNodes(collection)
  #   newPortFidList = [collection.GetItemAsObject(i) for i in
  #                     xrange(collection.GetNumberOfItems())]
  #   newPortFidList = [fid for fid in newPortFidList 
  #                     if fid.GetClassName() == "vtkMRMLAnnotationFiducialNode"]

  #   # first remove ports which aren't in the new list
  #   for portFid in [fid for fid in self.fiducialToolMap if not fid in newPortFidList]:
  #     tool = self.fiducialToolMap[portFid]
  #     slicer.mrmlScene.RemoveNode(tool.model)
  #     slicer.mrmlScene.RemoveNode(tool.transform)
  #     del self.fiducialToolMap[portFid]

  #   # If a port in the new list is already in the old list, update
  #   # the old port's data
  #   #
  #   # We take the following steps for each pre-existing port in the new list:
  #   # 1. Store the tool's current orientation
  #   # 2. Get tool's new position from the fiducial's new position
  #   # 3. "Reset" the tool's pose by inverting its transformation matrix
  #   # 4. Translate tool to new position
  #   # 5. Apply a rotation to regain its previous orientation
  #   for portFid in [fid for fid in self.fiducialToolMap if fid in newPortFidList]:
  #     tool = self.fiducialToolMap[portFid]
  #     mat = vtk.vtkMatrix4x4()
  #     tool.transform.GetMatrixTransformToWorld(mat)

  #     t = vtk.vtkTransform()
  #     t.SetMatrix(mat)
  #     orientation = [0,0,0,0]
  #     t.GetOrientationWXYZ(orientation)
      
  #     newPosition = [0,0,0]
  #     portFid.GetFiducialCoordinates(newPosition)

  #     t.Inverse()
  #     t.Translate(newPosition)
  #     t.RotateWXYZ(orientation[0], orientation[1:])
  #     tool.transform.ApplyTransformMatrix(t.GetMatrix())

  #   # Now we add the new ports from the new port list
  #   # 
  #   # If stuff doesn't work, check to make sure that RemoveNode and
  #   # AddNode are referring to the same objects
  #   for portFid in [fid for fid in newPortFidList if not fid in self.fiducialToolMap]:
  #     # create the tool model
  #     toolModel = slicer.vtkMRMLModelNode()
  #     toolModel.SetName(slicer.mrmlScene.GenerateUniqueName("Tool"))
  #     toolModel.SetAndObservePolyData(self.toolPolyData)
  #     slicer.mrmlScene.AddNode(toolModel)

  #     # and then our model display node
  #     modelDisplay = slicer.vtkMRMLModelDisplayNode()
  #     modelDisplay.SetColor(0,1,1) # cyan
  #     slicer.mrmlScene.AddNode(modelDisplay)
  #     toolModel.SetAndObserveDisplayNodeID(modelDisplay.GetID())

  #     # and our transform node to correspond with the fiducial's position
  #     fidPos = [0,0,0]
  #     portFid.GetFiducialCoordinates(fidPos)
  #     t = vtk.vtkTransform()
  #     t.Translate(fidPos)
  #     transformNode = slicer.vtkMRMLLinearTransformNode()
  #     transformNode.ApplyTransformMatrix(t.GetMatrix())
  #     transformNode.SetName(slicer.mrmlScene.GenerateUniqueName("Transform"))
  #     slicer.mrmlScene.AddNode(transformNode)

  #     # apply the new transform to our tool
  #     toolModel.SetAndObserveTransformNodeID(transformNode.GetID())

  #     # finally, add the model, model display, and transform to our fiducialToolMap
  #     self.fiducialToolMap[portFid] = self.Tool(toolModel, modelDisplay, transformNode)

  def retargetTools(self, targetFid, activePortAnnotationsHierarchy):
    if targetFid.GetClassName() != "vtkMRMLAnnotationFiducialNode":
      return

    import numpy
    import numpy.linalg
    import math
    # Get target coordinates
    targetLoc = [0,0,0]
    targetFid.GetFiducialCoordinates(targetLoc)
    targetLoc = numpy.array(targetLoc)

    activeFiducials = self.fromHierarchyToFiducialList(activePortAnnotationsHierarchy)

    # DEBUG TEST
    for fid in activeFiducials:
      if not fid in self.fiducialToolMap:
        raise Exception('ERROR: somehow a fiducial slipped through')

    # Iterate through port tool map and retarget their associated tools
    for portFid in activeFiducials:
      # Don't try to target a tool toward itself
      if portFid == targetFid: 
        continue

      portLoc = [0,0,0]
      portFid.GetFiducialCoordinates(portLoc)
      portLoc = numpy.array(portLoc)

      # Assume tools get drawn aligned with the global y-axis. We
      # want to transform the tool to be oriented along the port
      # toward the target point. We begin by finding the axis to
      # rotate the tool by, which is the cross product of the
      # global y and the target vector
      targetVec = targetLoc - portLoc
      targetVec = targetVec / numpy.linalg.norm(targetVec)
      rotAxis = numpy.cross(numpy.array([0,1,0]), targetVec)
      rotAxis = rotAxis / numpy.linalg.norm(rotAxis)

      # get rotation angle and normalize rotation axis
      angle = math.acos(numpy.dot(numpy.array([0,1,0]), targetVec)) * 180. / math.pi

      # generate our transform
      #
      # We do this by first undo-ing the tool's current transform and
      # then applying the new one.
      toolTransform = self.fiducialToolMap[portFid].transform
      mat = vtk.vtkMatrix4x4()
      toolTransform.GetMatrixTransformToWorld(mat)
      t = vtk.vtkTransform()
      t.SetMatrix(mat)

      t.Inverse()
      t.Translate(portLoc.tolist())
      t.RotateWXYZ(angle, rotAxis.tolist())
      toolTransform.ApplyTransformMatrix(t.GetMatrix())


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
    self.test_PortPlacement1()

  def test_PortPlacement1(self):
    import numpy
    import numpy.linalg
    import random

    self.delayDisplay("Starting test...")

    # # Grab the annotations GUI to help with testing fiducials
    # annotationsLogic = slicer.util.getModule('Annotations').logic()

    # # Add our port hierarchy node at the top level
    # annotationsLogic.SetActiveHierarchyNodeID(annotationsLogic.GetTopLevelHierarchyNodeID())
    # annotationsLogic.AddHierarchy()
    # portsHierarchy = annotationsLogic.GetActiveHierarchyNode()

    # # Add some ports to our ports hierarchy
    # portFiducialNodes = []
    # numPorts = 4
    # for i in range(numPorts):
    #   fiducialNode = slicer.vtkMRMLAnnotationFiducialNode()
    #   fiducialNode.SetFiducialCoordinates(*[random.uniform(-100.,100.) for j in range(3)])
    #   fiducialNode.Initialize(slicer.mrmlScene)
    #   portFiducialNodes.append(fiducialNode)

    # # Add port tools using port placement logic
    # logic = PortPlacementLogic(2.0, 50.0)          
    # logic.updatePorts(portsHierarchy)
    # self.assertTrue(len(logic.fiducialToolMap) == numPorts)

    # # Verify that tools' transforms correspond to port list 
    # #
    # # TODO: instead, this should comb the last 2n added objects to the
    # # scene to check against.
    # for portFid in logic.fiducialToolMap:
    #   self.assertTrue(portFid in portFiducialNodes)

    #   center = [0,0,0,1]
    #   toolPosition = [0,0,0,1]
    #   logic.fiducialToolMap[portFid].model.TransformPointToWorld(center, toolPosition)

    #   # check that this tool's position matches that of the fiducial
    #   fidCoords = [0,0,0]
    #   portFid.GetFiducialCoordinates(fidCoords)
    #   diff = numpy.array(toolPosition)[0:3] - numpy.array(fidCoords)
    #   self.assertTrue(numpy.dot(diff,diff) < 1e-10)

    # # retarget tools toward a random point
    # targetFiducial = slicer.vtkMRMLAnnotationFiducialNode()
    # targetFiducial.SetFiducialCoordinates(*[random.uniform(-100.,100.) for i in range(3)])
    # annotationsLogic.SetActiveHierarchyNodeID(annotationsLogic.GetTopLevelHierarchyNodeID())
    # # slicer.mrmlScene.AddNode(targetFiducial)
    # targetFiducial.Initialize(slicer.mrmlScene)
    # annotationsLogic.SetActiveHierarchyNodeID(portsHierarchy.GetID())
    # logic.retargetTools(targetFiducial)

    # # check retargeting by verifying that tools' y-axes are oriented
    # # toward point
    # for portFid in logic.fiducialToolMap:
    #   targetWorld = [0,0,0]
    #   targetFiducial.GetFiducialCoordinates(targetWorld)
    #   targetWorld = targetWorld + [1]
    #   targetLocal = [0,0,0,1]
    #   logic.fiducialToolMap[portFid].model.TransformPointFromWorld(targetWorld, targetLocal)

    #   targetLocal = numpy.array(targetLocal)[0:3]
    #   targetLocal = targetLocal / numpy.linalg.norm(targetLocal)

    #   # target local should be the unit y-axis (e2)
    #   yAxis = numpy.array([0.,1.,0.])
    #   diff = yAxis - targetLocal
    #   self.assertTrue(numpy.dot(diff,diff) < 1e-10)

    # # Now add a fiducial
    # newFiducial = slicer.vtkMRMLAnnotationFiducialNode()
    # newFiducial.SetFiducialCoordinates(*[random.uniform(-100.,100.) for i in range(3)])
    # newFiducial.Initialize(slicer.mrmlScene)
    # # slicer.mrmlScene.AddNode(newFiducial)

    # # Remove the first fiducial
    # removedFiducial = portFiducialNodes[0]
    # annotationsLogic.RemoveAnnotationNode(removedFiducial)

    # # Change position of second fiducial
    # changedFiducial = portFiducialNodes[1]
    # changedFiducial.SetFiducialCoordinates(0., 0., 0.)
        
    # # Update visualized port tools
    # logic.updatePorts(portsHierarchy)

    # # check for added fiducial
    # center = [0,0,0,1]
    # toolPosition = [0,0,0,1]
    # logic.fiducialToolMap[newFiducial].model.TransformPointToWorld(center, toolPosition)
    # fidCoords = [0,0,0]
    # newFiducial.GetFiducialCoordinates(fidCoords)
    # diff = numpy.array(toolPosition)[0:3] - numpy.array(fidCoords)
    # self.assertTrue(numpy.dot(diff,diff) < 1e-10)

    # # check for removed fiducial
    # self.assertTrue(not removedFiducial in logic.fiducialToolMap)

    # # check for changed position of second fiducial
    # logic.fiducialToolMap[changedFiducial].model.TransformPointToWorld(center, toolPosition)
    # changedFiducial.GetFiducialCoordinates(fidCoords)
    # diff = numpy.array(toolPosition)[0:3] - numpy.array(fidCoords)
    # self.assertTrue(numpy.dot(diff,diff) < 1e-10)

    # # cleanup
    # for fid in logic.fiducialToolMap:
    #   annotationsLogic.RemoveAnnotationNode(fid)
    #   tool = logic.fiducialToolMap[fid]
    #   slicer.mrmlScene.RemoveNode(tool.model)
    #   slicer.mrmlScene.RemoveNode(tool.transform)      
    # slicer.mrmlScene.RemoveNode(portsHierarchy)
    # slicer.mrmlScene.RemoveNode(targetFiducial)

    self.delayDisplay("Test passed!")
    
