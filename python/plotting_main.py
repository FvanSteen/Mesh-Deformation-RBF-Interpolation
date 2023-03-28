# plotting_main

import os
from bdryScatter import bdryScatterFuns
from funs import getMeshQualParams, getPlotData
from meshQualPlot import meshQualPlot

meshPath = os.path.dirname(os.path.abspath(__file__)) + "/Meshes/"

#%% ----- 2D boundary scatter plot -----
if(0):
    fNameInit = "25x25.su2"
    fName = ["25x25_def.su2"]
    
    [f_init,v_init,elemType,bdryPnts_init,_,_,_] = getPlotData(meshPath + fNameInit)
    
    # Only the deformed mesh
    bdryScatterFuns.bdryScatter(meshPath + fName[0])
    
    # Deformed mesh and initial mesh
    bdryScatterFuns.bdryScatter2(v_init,bdryPnts_init, meshPath + fName[0])

#%% ----- 2D mesh quality plot -----

if(0):
    fNameInit = "25x25.su2"
    fName = ["25x25_def.su2"]
    
    [f_init,v_init,elemType,bdryPnts_init,_,_,_] = getPlotData(meshPath + fNameInit)
    
    [alphas_0,_,_,_] = getMeshQualParams(f_init,v_init,elemType)
    
    meshQualPlot(meshPath, fName,fNameInit,alphas_0)


#%% ----- 3D BOUNDARY SCATTER PLOT -----

if(1):
#    fName = "stator_per_ffd.su2"
    fName = "stator_per_ffd_deform.su2"
    
    [f_init,v_init,elemType,bdryPnts_init,markerTags, nElemsMarks, FFD_pnts] = getPlotData(meshPath + fName)
    
    print("Available marker tags for plotting: ", markerTags, "(excluding SEND_RECEIVE)")
    
    plotTag = ["BLADE", "HUB","INFLOW","OUTFLOW","PER1","PER2", "SHROUD"]   
    
    bdryScatterFuns.bdryScatter3D(v_init, bdryPnts_init, markerTags, nElemsMarks, plotTag, FFD_pnts)
    
