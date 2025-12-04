# trace generated using paraview version 6.0.0-RC1
#import paraview
#paraview.compatibility.major = 6
#paraview.compatibility.minor = 0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


import argparse

p = argparse.ArgumentParser(description="Extract fault variables from a xdmf output file")
p.add_argument("--path",help="Input data folder ")                 # positional
p.add_argument("--prefix", help="Input xdmf filename")
p.add_argument("--outputdir", help="Input xdmf filename")
p.add_argument("-v", "--verbose", action="store_true",        # boolean flag
                   help="Increase logging verbosity")
args = p.parse_args()

if args.verbose:
        print(f"[info] reading {args.path} {args.prefix} ")


folder=args.path
prefix=args.prefix
outputdir=args.outputdir

print(f"{folder}{prefix}-fault.xdmf")
# create a new 'XDMF Reader'
jp4bFfaultxdmf = XDMFReader(registrationName=f"{prefix}-fault.xdmf", FileNames=[f"{folder}{prefix}-fault.xdmf"])

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# show data in view
jp4bFfaultxdmfDisplay = Show(jp4bFfaultxdmf, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
jp4bFfaultxdmfDisplay.Representation = 'Surface'

# reset view to fit data
renderView1.ResetCamera(False, 0.9)

# get the material library
materialLibrary1 = GetMaterialLibrary()

# show color bar/color legend
jp4bFfaultxdmfDisplay.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color transfer function/color map for 'partition'
partitionLUT = GetColorTransferFunction('partition')

# get opacity transfer function/opacity map for 'partition'
partitionPWF = GetOpacityTransferFunction('partition')

# get 2D transfer function for 'partition'
partitionTF2D = GetTransferFunction2D('partition')

# create a new 'Calculator'
calculator1 = Calculator(registrationName='Calculator1', Input=jp4bFfaultxdmf)

# Properties modified on calculator1
calculator1.Set(
    ResultArrayName='rake',
    Function='asin(Sld/sqrt(Sld^2+Sls^2))*180/3.1415926',
)

# show data in view
calculator1Display = Show(calculator1, renderView1, 'UnstructuredGridRepresentation')


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
Hide(jp4bFfaultxdmf, renderView1)

# show color bar/color legend
cellSize1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set scalar coloring
ColorBy(cellSize1Display, ('CELLS', 'Area'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(partitionLUT, renderView1)

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

# get color legend/bar for areaLUT in view renderView1
areaLUTColorBar = GetScalarBar(areaLUT, renderView1)


# trace defaults for the display properties.
calculator1Display.Representation = 'Surface'

# hide data in view
Hide(jp4bFfaultxdmf, renderView1)

# show color bar/color legend
calculator1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# get color transfer function/color map for 'rake'
rakeLUT = GetColorTransferFunction('rake')

# get opacity transfer function/opacity map for 'rake'
rakePWF = GetOpacityTransferFunction('rake')

# get 2D transfer function for 'rake'
rakeTF2D = GetTransferFunction2D('rake')

# Define variables to extract
variables_to_extract = ['fault-tag','Pn0',  'Td0', 'Ts0' ]

# Extract data at TimeStep = 0 (initial state)
if args.verbose:
    print("[info] Extracting data at TimeStep=0 (initial state)")

animationScene1.GoToFirst()

# rescale color and/or opacity maps used to exactly fit the current data range
calculator1Display.RescaleTransferFunctionToDataRange(False, True)

# Update view to ensure we're at the first time step
renderView1.Update()

# Save data at TimeStep=0
SaveData(f"{outputdir}stress_{prefix}_t0.csv", proxy=calculator1, ChooseArraysToWrite=1,
    CellDataArrays=variables_to_extract,
    FieldAssociation='Cell Data',
    AddTimeStep=1)

if args.verbose:
    print(f"[info] Saved initial state to stress_{prefix}_t0.csv")

# Extract data at last TimeStep (final state)
if args.verbose:
    print("[info] Extracting data at last TimeStep (final state)")

animationScene1.GoToLast()

# rescale color and/or opacity maps used to exactly fit the current data range
calculator1Display.RescaleTransferFunctionToDataRange(False, True)

# Update view to ensure we're at the last time step
renderView1.Update()

variables_to_extract = ['fault-tag','RT','PSR','ASl', 'SRs','SRd', 'Pn0', 'Sld', 'Sls', 'Td0', 'Ts0', 'Vr', 'rake','Area']

# Save data at last TimeStep
SaveData(f"{outputdir}stress_{prefix}_final.csv", proxy=cellSize1, ChooseArraysToWrite=1,
    CellDataArrays=variables_to_extract,
    FieldAssociation='Cell Data',
    AddTimeStep=1)

if args.verbose:
    print(f"[info] Saved final state to stress_{prefix}_final.csv")

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

# get layout
layout1 = GetLayout()

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(2190, 1382)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.Set(
    CameraPosition=[1769349.7134823536, 5051315.480176339, 360165.1414230793],
    CameraFocalPoint=[1790665.875, 5417453.999999997, -15000.000000000033],
    CameraViewUp=[-0.010655169861531415, 0.7159258631863311, 0.6980949976730451],
    CameraParallelScale=135790.16562242503,
)


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
