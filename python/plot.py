import os
from bdryScatter import bdryScatterFuns
from meshQualPlot import meshQualPlot
from funs import getMeshQualParams,getPlotData, getMeshQualParams3D, getMeshQuals3D
import matplotlib.pyplot as plt
import numpy as np

import matplotlib.collections
#plt.close('all')
# Setting directory to find .su2 files
figPath = os.path.dirname(os.path.abspath(__file__)) + "/figs/"
os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\Meshes')

if(0):
    fNameInit = '/25x25.su2'
    #fNameInit = '/turbine_row.su2'
    #fNameInit = '/mesh_NACA0012_inv.su2'
    #fNameInit = '/gis_nthx001_mesh.su2'
    #fileNames = ['/25x25_def_none.su2','/25x25_def_periodic.su2','/25x25_def_fixed.su2','/25x25_def_moving.su2']
    fileNames = ['/25x25_def.su2']
    
    #fileNames = ['/mesh_NACA0012_inv_def.su2']
    [f_init,v_init,elemType,bdryPnts_init,_,_,_] = getPlotData(fNameInit)
    
    [alphas_0,_,_,_] = getMeshQualParams(f_init,v_init,elemType)
#%%
#a = v_init[98320]
#b = v_init[97661]
#qc = a +0.25*(b-a)
#%%

#bdryPnts = np.unique(bdryPnts_init)
#blockidx = np.array([49,50,51,52,53,54,57,58,59,60,61,62])
#periodicIdx = np.append(np.arange(1,25), np.arange(111-11-25,111-12))
#blockPnts = bdryPnts[blockidx]
#
#bdryPnts = np.delete(bdryPnts, blockidx)
#periodicPnts = bdryPnts[periodicIdx]
#bdryPnts = np.delete(bdryPnts,periodicIdx)
#
#pc = matplotlib.collections.PolyCollection(v_init[f_init],cmap='seismic', facecolors="white", edgecolor="black",linewidth=0.4)
#
#labels= ["_nolegend_", "control", "sliding", "periodic"]
##labels= [ "control", "internal"]
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#polys = ax.add_collection(pc)
##plt.plot([0,1,1,0,0],[0,0,1,1,0], color = "black", linewidth=0.5)
#plt.scatter(v_init[blockPnts][:,0],v_init[blockPnts][:,1], s = 15)
#plt.scatter(v_init[bdryPnts][:,0],v_init[bdryPnts][:,1], s = 15)
#plt.scatter(v_init[periodicPnts][:,0],v_init[periodicPnts][:,1], s = 15)
#
#
#ax.set_xlim([-0.05,1.05])
#ax.set_ylim([-0.05,1.05])
#plt.axis('equal')
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.legend(labels, bbox_to_anchor=(1.04,1), loc='upper left')
#plt.show()
#plt.savefig(figPath + "per_3.png", dpi=800,bbox_inches='tight')
##%%
#
#bdryPnts = np.unique(bdryPnts_init)
#internalPnts = np.arange(0,np.size(v_init,0))
#internalPnts = np.delete(internalPnts,bdryPnts)
#
#
#pc = matplotlib.collections.PolyCollection(v_init[f_init],cmap='seismic', facecolors="white", edgecolor="black",linewidth=0.5)
#
#labels= ["_nolegend_", "control", "internal"]
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#polys = ax.add_collection(pc)
#plt.scatter(v_init[bdryPnts][:,0],v_init[bdryPnts][:,1], s = 15)
#plt.scatter(v_init[internalPnts][:,0],v_init[internalPnts][:,1], s= 15)
#
#ax.set_xlim([-0.05,1.05])
#ax.set_ylim([-0.05,1.05])
#plt.axis('equal')
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.legend(labels, bbox_to_anchor=(1.04,1), loc='upper left')
#plt.show()
#
#
#plt.savefig(figPath + "rbf_def.png", dpi=800,bbox_inches='tight')


#[f,v,elemType,_] = getPlotData('/25x25.su2')  
#[f_init, v_init, elemType, bdryPnts_init] = getPlotData('/mesh_NACA0012_inv.su2')
#[f_init,v_init,elemType,bdryPnts_init] = getPlotData('/turbine_row.su2')
#[f,v,elemType] = getPlotData('/su2mesh.su2')
#[f,v,elemType] = getPlotData('/mesh_original.su2')




# Provide filenames of the meshes
#fileNames = ['/25x25_def_ps_none_noCorrect_DataRed_e1.su2','/25x25_def_ps_none_noCorrect_DataRed_e2.su2','/25x25_def_ps_none_noCorrect_DataRed_e3.su2','/25x25_def_ps_none_noCorrect_DataRed_e4.su2','/25x25_def_ps_none_correct_DataRed_e1.su2','/25x25_def_ps_none_correct_DataRed_e2.su2','/25x25_def_ps_none_correct_DataRed_e3.su2','/25x25_def_ps_none_correct_DataRed_e4.su2']
#fileNames = ['/375x375_def.su2']
#fileNames = ['/mesh_NACA0012_inv_def.su2']
#fileNames = ['/LS89_turbine_def.su2']
#fileNames= ['/25x25_def_ref.su2','/25x25_def_err1.su2','/25x25_def_err2.su2','/25x25_def_err3.su2','/25x25_def_err4.su2','/25x25_def_err5.su2']
#fileNames = ['/mesh_NACA0012_inv_def_none_ref.su2','/mesh_NACA0012_inv_def_none_greedy.su2','/mesh_NACA0012_inv_def_ps_ref.su2','/mesh_NACA0012_inv_def_ps_greedy.su2']
#fileNames = ['/su2mesh_def_ps_ref.su2']
#fileNames = ['/mesh_original.su2']
#fileNames = ['/turbine_row_def.su2']


#fileNames = ['/gis_nthx001_mesh_def.su2']
#graphNames = ['non-periodic', 'periodic in y', 'periodic in y,\nmoving boundaries, fixed vertices', 'periodic in y,\nmoving boundaries, moving vertices']
#graphNames = ['','Regular RBF', 'Greedy 1e-1', 'Greedy 1e-2', 'Greedy 1e-3', 'Greedy 1e-4', 'Greedy 1e-5']
#graphNames = ['standard','periodic', 'fixed', 'moving']


#%% Mesh Quality plots
    #import numpy as np
    graphNames = ['','']
    #fileNames = ['/turbine_row_d.su2']
    meshQualPlot(fileNames,fNameInit,graphNames,alphas_0)
    #plt.scatter(qc[0],qc[1])
    #x = np.argwhere(np.isnan(meshQual))
    #plt.savefig(figPath + "turbine_row_ps_fixed_zoom.png", dpi=800,bbox_inches='tight')
#%% MAKING A SCATTER PLOT OF THE BOUNDARY POINTS

### with the initial mesh ###
if(0):
#    bdryScatterFuns.bdryScatter2(v_init,bdryPnts_init, fileNames[0])

### without the initial mesh ###
    bdryScatterFuns.bdryScatter(fileNames[0])



#%%
if(1):
    #import numpy as np
    [f_init,v_init,elemType,bdryPnts_init,markerTags, nElemsMarks, FFD_pnts] = getPlotData("/9x9x9_def.su2")
    plotTag= ["FRONT","BACK","LOWER","UPPER","LEFT","RIGHT", "BLOCK"]

    bdryScatterFuns.bdryScatter3D(v_init, bdryPnts_init, markerTags, nElemsMarks, plotTag, FFD_pnts)


#%%
#pnts = []
#for i in range(np.shape(v_init)[0]): 
#    if v_init[i,0] < 1.5 and v_init[i,0] > -0.5 and v_init[i,2] < 0.25 and v_init[i,2] > -0.25 and v_init[i,1] >= 0 and v_init[i,1] <1.3:
#        pnts.append(i)
#    
#fig = plt.figure()
#ax = fig.add_subplot(projection='3d')
#ax.scatter(v_init[pnts,0],v_init[pnts,1],v_init[pnts,2])
# 
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.set_zlabel('z')
#%%
#from mpl_toolkits import mplot3d
#
#fig = plt.figure()
#ax = plt.axes(projection="3d")

#ax.scatte3D()
#fname = '/25x25x5_def.su2'
#bdryScatterFuns.bdryScatter3D(fname)
#%%
if(1):
    import numpy as np
    import math
    import matplotlib.collections
    from colMap import colMap
    
    cmapMatlab = colMap()
    fileNames = ['/9x9x9_def.su2']
    
    [f_init,v_init,elemType,bdryPnts_init,_,_,_] = getPlotData('/9x9x9.su2')
    
    
    graphNames = ["",""]
    cutAxis = 2
    cutPlaneLoc = 0.5
    
    dims = [0,1,2]
    dims.remove(cutAxis)
    [f,v,elemType,_, markerTags, nElemsMarks,FFD_pnts] = getPlotData(fileNames[0]) 
    
    idxCutElems = np.array([],dtype = 'int')
    bdryPnts_init_unique = np.unique(bdryPnts_init)
    for i in range(len(f)):
        if np.any(v[f][i][:,cutAxis] < cutPlaneLoc) and np.any(v[f][i][:,cutAxis] >= cutPlaneLoc):
            if set(f[i]).issubset(bdryPnts_init_unique) == False:
                idxCutElems = np.append(idxCutElems, i)
            7
    
    [alphas_0,_,_,_] = getMeshQualParams3D(f_init[idxCutElems,:],v_init,elemType[idxCutElems])

    #[f,v,elemType,_,_,_,_] = getPlotData(fileNames[0])
    meshQual = getMeshQuals3D(f[idxCutElems],v,alphas_0,elemType[idxCutElems])
    print("Min mesh quality: \t", round(np.min(meshQual),5)) 
    print("Max mesh quality: \t", round(np.max(meshQual),5)) 
    print("Mean mesh quality: \t", round(np.mean(meshQual),5))
    
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
    
        midpoint = np.sum(cutLoc[0:cnt,:],axis=0)/cnt
    
        angles = np.empty(cnt,dtype=float)
        for ii in range(cnt):
            angles[ii] = math.atan2(cutLoc[ii,1]-midpoint[1],cutLoc[ii,0]-midpoint[0])*180/np.pi
    
        cutLoc[0:cnt,:] = cutLoc[0:cnt,:][angles.argsort()]
    
        v_cut[i][0:cnt,:] = cutLoc[0:cnt,:]
    #    print(cutLoc[0:cnt,:])
    
    #    plt.scatter(cutLoc[0:cnt,0],cutLoc[0:cnt,1], color= "blue") 
                
    #    for ii in range(len(v1)):
                
    #        print(v[f][idxCutElems[i]][v1,cutAxis], v[f][idxCutElems[i]][v2,cutAxis])
            
    colors = cmapMatlab(plt.Normalize(0,1)(meshQual))
    fig = plt.figure()  
    pc = matplotlib.collections.PolyCollection(v_cut,cmap=cmapMatlab, facecolors=colors, edgecolor="black",linewidth=0.1)
    ax = fig.add_subplot(1,1,1)
    
    polys = ax.add_collection(pc)
    pc.set_array(None)
    ax.autoscale()  
    ax.set_aspect('equal')
    polys.set_clim(0,1)
    plt.colorbar(polys, ax=ax, shrink=1.0/(np.size(graphNames)-1))
    plt.show()
#%%
if(0):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    cutElems = v[f[idxCutElems]]
    for i in range(len(cutElems)):
        ax.scatter(cutElems[i][:,0],cutElems[i][:,1],cutElems[i][:,2],marker = "1", color='blue')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_zlim(0,1)
    #plotTag = ["BLOCK", "LOWER","UPPER","FRONT","BACK","LEFT","RIGHT"]
    #bdryScatterFuns.bdryScatter3D(v, bdryPnts_init, markerTags, nElemsMarks, plotTag, FFD_pnts)
    

#%%
#import os 
#import numpy as np
#os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\Classical RBF\\defs')
#x = os.getcwd()
#o = open(x + "\\gis_nthx001_mesh_deformation.txt", "w")
#i_base = open(x + "\\MoveSurface_baseline.txt", "r") 
#i_def = open(x + "\\MoveSurface_deformed.txt", "r") 
#
#lines_base = i_base.read().splitlines()
#lines_def = i_def.read().splitlines()
#
#delta = np.empty((len(lines_base), 2), dtype = float)
#index = np.empty((len(lines_base),1), dtype = int)
#for i in range(0,len(lines_base)):
#    i_split_base = lines_base[i].split()
#    i_split_def = lines_def[i].split()
#    delta[i,0] = float(i_split_def[1]) - float(i_split_base[1])
#    delta[i,1] = float(i_split_def[2]) - float(i_split_base[2])
#    index[i] = int(i_split_def[0])
##    print(i_split_base)
##    delta[i,0] = i_split[]
##    print(i.split()[0])
#
#
#for idx, c1,c2 in zip(index[:,0],delta[:,0],delta[:,1]):
#    o.write("{0}\t{1}\t{2}".format(idx, c1,c2) + "\n")
#o.close()

