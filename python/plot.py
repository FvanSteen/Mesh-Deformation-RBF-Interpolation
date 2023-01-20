import os
from bdryScatter import bdryScatterFuns
from meshQualPlot import meshQualPlot
from funs import getMeshQualParams,getPlotData
import matplotlib.pyplot as plt

#plt.close('all')
# Setting directory to find .su2 files
os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\Classical RBF\\Meshes')

#fNameInit = '/LS89_turbine.su2'
#fNameInit = '/mesh_NACA0012_inv.su2'
#fNameInit = '/25x25.su2'
#fNameInit = '/turbine_row.su2'
fNameInit = '/5x5x5.su2'

[f_init,v_init,elemType,bdryPnts_init] = getPlotData(fNameInit)

#[f,v,elemType,_] = getPlotData('/25x25.su2')  
#[f_init, v_init, elemType, bdryPnts_init] = getPlotData('/mesh_NACA0012_inv.su2')
#[f_init,v_init,elemType,bdryPnts_init] = getPlotData('/turbine_row.su2')
#[f,v,elemType] = getPlotData('/su2mesh.su2')
#[f,v,elemType] = getPlotData('/mesh_original.su2')
[alphas_0,_,_,_] = getMeshQualParams(f_init,v_init,elemType)





# Provide filenames of the meshes
#fileNames = ['/25x25_def_ps_none_noCorrect_DataRed_e1.su2','/25x25_def_ps_none_noCorrect_DataRed_e2.su2','/25x25_def_ps_none_noCorrect_DataRed_e3.su2','/25x25_def_ps_none_noCorrect_DataRed_e4.su2','/25x25_def_ps_none_correct_DataRed_e1.su2','/25x25_def_ps_none_correct_DataRed_e2.su2','/25x25_def_ps_none_correct_DataRed_e3.su2','/25x25_def_ps_none_correct_DataRed_e4.su2']
#fileNames = ['/25x25_def.su2']
#fileNames = ['/mesh_NACA0012_inv_def.su2']
#fileNames = ['/LS89_turbine_def.su2']
#fileNames= ['/25x25_def_ref.su2','/25x25_def_err1.su2','/25x25_def_err2.su2','/25x25_def_err3.su2','/25x25_def_err4.su2','/25x25_def_err5.su2']
#fileNames = ['/mesh_NACA0012_inv_def_none_ref.su2','/mesh_NACA0012_inv_def_none_greedy.su2','/mesh_NACA0012_inv_def_ps_ref.su2','/mesh_NACA0012_inv_def_ps_greedy.su2']
#fileNames = ['/su2mesh_def_ps_ref.su2']
#fileNames = ['/mesh_original.su2']
#fileNames = ['/turbine_row_def.su2']
#fileNames = ['/25x25x25.su2']

#graphNames = ['non-periodic', 'periodic in y', 'periodic in y,\nmoving boundaries, fixed vertices', 'periodic in y,\nmoving boundaries, moving vertices']
#graphNames = ['','Regular RBF', 'Greedy 1e-1', 'Greedy 1e-2', 'Greedy 1e-3', 'Greedy 1e-4', 'Greedy 1e-5']
graphNames = ['','','','']


#%% Mesh Quality plots
#import numpy as np
#meshQualPlot(fileNames,fNameInit,graphNames,alphas_0)

#x = np.argwhere(np.isnan(meshQual))

#%% MAKING A SCATTER PLOT OF THE BOUNDARY POINTS

### with the initial mesh ###
#bdryScatterFuns.bdryScatter2(v_init,bdryPnts_init, fileNames[0])

### without the initial mesh ###
#bdryScatterFuns.bdryScatter(fileNames[0])

#bdryScatterFuns.bdryScatter3D(fileNames[0])

#%%
#from mpl_toolkits import mplot3d
#
#fig = plt.figure()
#ax = plt.axes(projection="3d")
#ax.scatte3D()
fname = '/5x5x5_def.su2'
bdryScatterFuns.bdryScatter3D(fname)
#%%
import numpy as np
import math
import matplotlib.collections
fileNames = ['/5x5x5_def.su2']
cutAxis = 2
cutPlaneLoc = 0.5

dims = [0,1,2]
dims.remove(cutAxis)
[f,v,elemType,_] = getPlotData(fileNames[0]) 

idxCutElems = np.array([],dtype = 'int')

for i in range(len(f)):
    if np.any(v[f][i][:,cutAxis] < cutPlaneLoc) and np.any(v[f][i][:,cutAxis] >= cutPlaneLoc):
        
        idxCutElems = np.append(idxCutElems, i)
        

v1 = [0,1,2,3,0,1,2,3,4,5,6,7]
v2 = [1,2,3,0,4,5,6,7,5,6,7,4]

verticesInd = np.array([[0,1],[1,2],[2,3],[3,0],[0,4],[1,5],[2,6],[3,7],[4,5],[5,6],[6,7],[7,4]])
v_cut = np.empty((len(idxCutElems),8,2),dtype = float)
v_cut[:] = np.nan
#plt.figure();
for i in range(len(idxCutElems)):
    edges = v[f][idxCutElems[i]][verticesInd]    
    cutLoc = np.empty((12,2), dtype = float)
    cnt = 0
    for ii in range(len(edges)):
        if np.any(edges[ii][:,cutAxis] < cutPlaneLoc) and np.any(edges[ii][:,cutAxis] >= cutPlaneLoc):
            delta = edges[ii][1,:] - edges[ii][0,:]
            cutLoc[cnt,:] = edges[ii][0,dims] + delta[dims]*((cutPlaneLoc-edges[ii][0,cutAxis])/delta[cutAxis])
            cnt = cnt+1

    midpoint = np.sum(cutLoc,axis=0)/np.size(cutLoc,axis=0)
    angles = np.empty(cnt,dtype=float)
    for ii in range(cnt):
#        print(ii) 
        angles[ii] = math.atan2(cutLoc[ii,1]-midpoint[1],cutLoc[ii,0]-midpoint[0])*180/np.pi
     
    cutLoc[0:cnt,:] = cutLoc[0:cnt,:][angles.argsort()]
    
    v_cut[i][0:cnt,:] = cutLoc[0:cnt,:]
#    print(cutLoc[0:cnt,:])

#    plt.scatter(cutLoc[0:cnt,0],cutLoc[0:cnt,1], color= "blue")
            
#    for ii in range(len(v1)):
            
#        print(v[f][idxCutElems[i]][v1,cutAxis], v[f][idxCutElems[i]][v2,cutAxis])
        

fig = plt.figure()  
pc = matplotlib.collections.PolyCollection(v_cut,cmap='seismic', facecolors='blue', edgecolor="black",linewidth=0.5)
ax = fig.add_subplot(1,1,1)

polys = ax.add_collection(pc)

