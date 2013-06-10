import os
import unittest
from __main__ import vtk, qt, ctk, slicer

#
# PortPlacement
#

class PortPlacement:
  def __init__(self, parent):
    parent.title = "Port Placement" # TODO make this more human readable by adding spaces
    parent.categories = ["Work in Progress"]
    parent.dependencies = []
    parent.contributors = ["Luis G. Torres (UNC)"] # replace with "Firstname Lastname (Org)"
    parent.helpText = """
    This is an example of scripted loadable module bundled in an extension.
    """
    parent.acknowledgementText = """
    This file was originally developed by Luis G. Torres (UNC).
""" # replace with organization, grant and thanks.
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
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

    #
    # radius spin box
    # 
    self.radiusSpinBox = qt.QDoubleSpinBox()
    self.radiusSpinBox.setMinimum(0.0)
    self.radiusSpinBox.setMaximum(10.0)
    self.radiusSpinBox.setValue(2.0)
    parametersFormLayout.addRow("Tool radius: ", self.radiusSpinBox)

    #
    # length spin box
    #
    self.lengthSpinBox = qt.QDoubleSpinBox()
    self.lengthSpinBox.setMinimum(0.0)
    self.lengthSpinBox.setMaximum(250.0)
    self.lengthSpinBox.setValue(150.0)
    parametersFormLayout.addRow("Tool length: ", self.lengthSpinBox)    

    #
    # target selector
    #
    self.targetSelector = slicer.qMRMLNodeComboBox()
    self.targetSelector.nodeTypes = ["vtkMRMLAnnotationFiducialNode"]
    self.targetSelector.addEnabled = False
    self.targetSelector.removeEnabled = False
    self.targetSelector.noneEnabled = False
    self.targetSelector.setMRMLScene(slicer.mrmlScene)
    self.targetSelector.setToolTip("Pick the surgical targets that the tools should try to reach.")
    parametersFormLayout.addRow("Target Point: ", self.targetSelector)

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
    parametersFormLayout.addRow("List of Ports: ", self.portListSelector)

    #
    # Update button
    #
    self.updateButton = qt.QPushButton("Update")
    self.updateButton.toolTip = "Update the port placement visualization"
    self.updateButton.enabled = True
    parametersFormLayout.addRow(self.updateButton)

    # connections
    self.updateButton.connect('clicked(bool)', self.onUpdateButton)
    self.radiusSpinBox.connect('valueChanged(double)', self.onToolShapeChanged)
    self.lengthSpinBox.connect('valueChanged(double)', self.onToolShapeChanged)

    # Add vertical spacer
    self.layout.addStretch(1)

    # set default values for tool shape
    self.onToolShapeChanged()

  def onUpdateButton(self):
    logic = PortPlacementLogic()
    logic.run(self.targetSelector.currentNode(), 
              self.portListSelector.currentNode(),
              self.toolPolyData)

  def onToolShapeChanged(self):
    toolSrc = vtk.vtkCylinderSource()
    toolSrc.SetRadius(self.radiusSpinBox.value)
    toolSrc.SetHeight(self.lengthSpinBox.value)
    self.toolPolyData = toolSrc.GetOutput()
    self.onUpdateButton()

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
  def __init__(self):
    pass

  def run(self,targetFid,portFidList,toolPolyData):
    import numpy
    import numpy.linalg
    import math

    if targetFid and portFidList:
      # Get target coordinates
      targetLoc = [0,0,0]
      targetFid.GetFiducialCoordinates(targetLoc)
      targetLoc = numpy.array(targetLoc)
      
      # Get port list
      portLocations = []
      collection = vtk.vtkCollection()
      portFidList.GetChildrenDisplayableNodes(collection)
      for fid in [collection.GetItemAsObject(i) for i in 
                  xrange(collection.GetNumberOfItems())]:
        coords = [0,0,0]
        fid.GetFiducialCoordinates(coords)
        portLocations.append(numpy.array(coords))

      # Let's make some cylinders!
      # We start with generating the polygonal data
      # toolSrc = vtk.vtkCylinderSource()
      # toolSrc.SetHeight(150)
      # toolSrc.SetRadius(2.0)
      # toolSrc.SetResolution(100)
      # toolPolyData = toolSrc.GetOutput()

      # Iterate through port locations and create the associated tools
      for portLoc in portLocations:
        # First add a tool model to our scene
        toolModel = slicer.vtkMRMLModelNode()
        toolModel.SetName(slicer.mrmlScene.GenerateUniqueName("Tool"))
        toolModel.SetAndObservePolyData(toolPolyData)
        slicer.mrmlScene.AddNode(toolModel)

        # and then our model display node
        modelDisplay = slicer.vtkMRMLModelDisplayNode()
        modelDisplay.SetColor(0,1,1) # cyan
        slicer.mrmlScene.AddNode(modelDisplay)
        toolModel.SetAndObserveDisplayNodeID(modelDisplay.GetID())
        
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
        t = vtk.vtkTransform()
        t.Translate(portLoc.tolist())
        t.RotateWXYZ(angle, rotAxis.tolist())
        transformNode = slicer.vtkMRMLLinearTransformNode()
        transformNode.ApplyTransformMatrix(t.GetMatrix())
        transformNode.SetName(slicer.mrmlScene.GenerateUniqueName("Transform"))
        slicer.mrmlScene.AddNode(transformNode)

        # apply the new transform to our tool
        toolModel.SetAndObserveTransformNodeID(transformNode.GetID())

    return True


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
    """ Ideally you should have several levels of tests.  At the lowest level
    tests sould exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    import urllib
    downloads = (
        ('http://slicer.kitware.com/midas3/download?items=5767', 'FA.nrrd', slicer.util.loadVolume),
        )

    for url,name,loader in downloads:
      filePath = slicer.app.temporaryPath + '/' + name
      if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
        print('Requesting download %s from %s...\n' % (name, url))
        urllib.urlretrieve(url, filePath)
      if loader:
        print('Loading %s...\n' % (name,))
        loader(filePath)
    self.delayDisplay('Finished with download and loading\n')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = PortPlacementLogic()
    self.assertTrue( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
