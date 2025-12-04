# trace generated using paraview version 5.13.3
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 13

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
welF400tapfaultxdmf = FindSource('welF400tap-fault.xdmf')

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=welF400tapfaultxdmf)

# Properties modified on calculator1
calculator1.Function = 'sqrt(Sls^2)'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
calculator1Display = Show(calculator1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
calculator1Display.Representation = 'Surface'

# hide data in view
Hide(welF400tapfaultxdmf, renderView1)

# show color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color transfer function/color map for 'Result'
resultLUT = GetColorTransferFunction('Result')

# get opacity transfer function/opacity map for 'Result'
resultPWF = GetOpacityTransferFunction('Result')

# get 2D transfer function for 'Result'
resultTF2D = GetTransferFunction2D('Result')

# Properties modified on calculator1
calculator1.Function = 'sqrt(Sls^2)\t'

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Cell Size'
cellSize1 = CellSize(registrationName='CellSize1', Input=calculator1)

# Properties modified on cellSize1
cellSize1.ComputeVertexCount = 0
cellSize1.ComputeLength = 0
cellSize1.ComputeVolume = 0

# show data in view
cellSize1Display = Show(cellSize1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
cellSize1Display.Representation = 'Surface'

# hide data in view
Hide(calculator1, renderView1)

# show color bar/color legend
cellSize1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(cellSize1Display, ('CELLS', 'Area'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(resultLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
cellSize1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
cellSize1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'Area'
areaLUT = GetColorTransferFunction('Area')

# get opacity transfer function/opacity map for 'Area'
areaPWF = GetOpacityTransferFunction('Area')

# get 2D transfer function for 'Area'
areaTF2D = GetTransferFunction2D('Area')

# save data
SaveData('/Users/DuoL/Documents/NSHM/Central/cov/test1.csv', proxy=cellSize1, ChooseArraysToWrite=1,
    CellDataArrays=['ASl', 'Area', 'Result'],
    FieldAssociation='Cell Data')

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

# get layout
layout1 = GetLayout()

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(2110, 1382)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [1773142.0381676953, 5259427.992119528, 4210.473945453535]
renderView1.CameraFocalPoint = [1761299.0466482756, 5429951.330107098, -14599.999999999674]
renderView1.CameraViewUp = [0.11993860171508974, 0.11707959411839043, 0.9858534883336885]
renderView1.CameraParallelScale = 44508.07141042609


##--------------------------------------------
## You may need to add some code at the end of this python script depending on your usage, eg:
#
## Render all views to see them appears
# RenderAllViews()
#
## Interact with the view, usefull when running from pvpython
# Interact()
#
## Save a screenshot of the active view
# SaveScreenshot("path/to/screenshot.png")
#
## Save a screenshot of a layout (multiple splitted view)
# SaveScreenshot("path/to/screenshot.png", GetLayout())
#
## Save all "Extractors" from the pipeline browser
# SaveExtracts()
#
## Save a animation of the current active view
# SaveAnimation()
#
## Please refer to the documentation of paraview.simple
## https://www.paraview.org/paraview-docs/latest/python/paraview.simple.html
##--------------------------------------------