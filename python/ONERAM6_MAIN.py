import os
from ONERAM6_funs import rewriteSU2File, getFFDdef, getOutputSU2File, getDefFile, scatteredWing, getDisp
from funs import getPlotData
from bdryScatter import bdryScatterFuns
import numpy as np
# CHANGING THE DIRECTORY
os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\Meshes')

# Editing the file for right format
# ONLY HAS TO BE DONE ONCE, PERFORM MANUAL CHECK TO SEE WHETHER ALL LINES HAVE BEEN DONE
fName = "\\mesh_ONERAM6_inv_FFD.su2"
#rewriteSU2File(fName)

# FILE DEDICATED TO ALL THE ADAPTATION FOR THE ONERA M6 WING

#%% OBTAINING THE DATA FROM THE INITIAL MESH FILE
fileNames = ['/OneraM6_fixed.su2']
[f_init,v_init,elemType,bdryPnts_init, markerTags, nElemsMarks,FFD_pnts] = getPlotData(fileNames[0])

#%% SCATTER PLOT OF THE DIFFERENT PARTS OF THE WING AND FFD BOX
plotTag= ["LOWER_SIDE","UPPER_SIDE", "TIP", "FFD"]

#fileNames = ['/OneraM6_fixed_wingonly_def.su2']
#[_,v_init,_,_,_,_,_] = getPlotData(fileNames[0])
#plotTag= ["XNORMAL_FACES", "YNORMAL_FACES","ZNORMAL_FACES","SYMMETRY_FACE"]
bdryScatterFuns.bdryScatter3D(v_init, bdryPnts_init, markerTags, nElemsMarks, plotTag, FFD_pnts)


#%% deformation of the FFD points 
# SINE FUNCTION FOR VERTICAL DISPLACEMENT FROM 0 TO MAC AT TIP
# SINE FUNCTION FOR WING TWIST FROM 0 AT ROOT TO 30 DEGS AT TIP
disp = getFFDdef(FFD_pnts)

# UNCOMMENT NEXT 2 LINES TO SEE THE DEFORMED FFD BOX
#newpoints = FFD_pnts + disp
#plotTag = ["LOWER_SIDE","UPPER_SIDE", "TIP", "FFD"]
#bdryScatterFuns.bdryScatter3D(v_init, bdryPnts_init, markerTags, nElemsMarks, plotTag, newpoints)

#%% Generating an SU2 file containing the wingpoints as internal points and the FFD box as boundary points which 
# can be given a prescribed displacement
wingTags = ["LOWER_SIDE","UPPER_SIDE", "TIP"]
getOutputSU2File(FFD_pnts, v_init, bdryPnts_init, markerTags, nElemsMarks, wingTags)

#%% GENERATING A DISPLACEMENT FILE FOR THE 
os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\defs')
defFileName = "\\ONERAM6_FFD_deformation.txt"
index = np.arange(0,np.shape(disp)[0])
getDefFile(defFileName, index, disp)

#%% AFTER DOING THE RBF INTERPOLATION WITH THE C++ CODE
os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\Meshes')

fileNames = ['/OneraM6_FFD.su2']
[_,v_FFD,_,bdryPnts_FFD, _, _,_] = getPlotData(fileNames[0])

fileNames = ['/OneraM6_FFD_def.su2']
[_,v_def,_,_, _, _,_] = getPlotData(fileNames[0])

scatteredWing(v_FFD,v_def, bdryPnts_FFD)

#%% 

os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\defs')
[index, disp] = getDisp(v_FFD, v_def, bdryPnts_FFD, bdryPnts_init, markerTags, wingTags, nElemsMarks)
getDefFile("\\ONERAM6_deformation.txt", index, disp)

#%% OBTAINING THE DATA FROM THE INITIAL MESH FILE
os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\Meshes')
fileNames = ['/OneraM6_fixed_def.su2']
[f_def,v_def,elemType,bdryPnts_def, markerTags, nElemsMarks,FFD_pnts] = getPlotData(fileNames[0])

# SCATTER PLOT OF THE DIFFERENT PARTS OF THE WING AND FFD BOX
plotTag= ["LOWER_SIDE","UPPER_SIDE", "TIP"]
bdryScatterFuns.bdryScatter3D(v_def, bdryPnts_def, markerTags, nElemsMarks, plotTag, FFD_pnts)

    