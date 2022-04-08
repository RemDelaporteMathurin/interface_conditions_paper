# trace generated using paraview version 5.7.0
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# get animation scene
animationScene1 = GetAnimationScene()

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

animationScene1.GoToNext()

# get active view
lineChartView1 = GetActiveViewOrCreate('XYChartView')
# uncomment following to set a specific view size
# lineChartView1.ViewSize = [1358, 928]

for i in range(0, 83):
    # export view
    ExportView('/home/rdelaporte/FESTIM/iter_1d_case/profiles_solute_c_continuity/{:.0f}.csv'.format(i), view=lineChartView1)

    animationScene1.GoToNext()
