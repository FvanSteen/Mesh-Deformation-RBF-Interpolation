import os
from bdryScatter import bdryScatterFuns
from meshQualPlot import meshQualPlot2D, meshQualPlot3D, nanElems
from funs import getMeshQualParams,getPlotData, getMeshQualParams3D, getMeshQuals3D, getMeshQual, getMeshQuals
import matplotlib.pyplot as plt
import numpy as np
from coordTransform import coordTransform


import matplotlib.collections
#plt.close('all')
# Setting directory to find .su2 files
figPath = os.path.dirname(os.path.abspath(__file__)) + "/figs/"
os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\Meshes')

#%%
#fileName = '/stator_per_ffd_mod_test_def.su2'
fileName = '/25x25_def.su2'
[f,v,elemType,bdryPnts,markerTags, nElemsMarks, FFD_pnts] = getPlotData(fileName)  

#v = coordTransform.toCylindrical(v)
#v = coordTransform.toPolar(v)
#%% 2D PLOTS

# Mesh quality
if(1):    
    graphName = ''    
    meshQualPlot2D(fileName, graphName, f, v, elemType)
#    plt.savefig(figPath + "25x25_init.png", dpi=800,bbox_inches='tight')
# Boundary scatter plots
if(0):
    # Without initial boundary
    init = 0
    bdryScatterFuns.bdryScatter(fileName, v, bdryPnts, init)
#    plt.savefig(figPath + "pseudo2.png", dpi=800,bbox_inches='tight')

#%%

#ax.legend(labels, bbox_to_anchor=(1.04,1), loc='upper left')
#plt.show()
#plt.savefig(figPath + "rbf_def.png", dpi=800,bbox_inches='tight')


#%% 3D PLOTS
    
    
# Mesh quality
if(0):
    graphName = ''
    cutAxis = 0
    cutPlaneLoc = 0.2995
    meshQualPlot3D(fileName, graphName, cutAxis, cutPlaneLoc, f, v, elemType)
#    plt.savefig(figPath + "stator_def2_0.1perc.png", dpi=800,bbox_inches='tight')
    
    
# Boundary scatter
if(0):

#    plotTag = ["BLADE","HUB","INFLOW","OUTFLOW","SHROUD"]
#    plotTag = ["HUB"]
    plotTag = []
#    plotTag = ["LEFT", "RIGHT", "UPPER","LOWER","FRONT","BACK","BLOCK"]
    bdryScatterFuns.bdryScatter3D(v, bdryPnts, markerTags, nElemsMarks, plotTag, FFD_pnts)
#    plt.savefig(figPath + "f1.png", dpi=800,bbox_inches='tight')

if(0):
    nanElems(fileName,v,f)
#    largerOneElems(fileName,v,f)
    
#%%
fName = "25x25"
data = np.genfromtxt('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\runHistory\\'+fName+ ".txt", skip_header=1, dtype=None, encoding='utf-8')

smode = 'none'
pmode = 'none'
greedy = 0
doubleEdge = 0

rows = ((data['f0']==smode) & (data['f1'] == pmode) & (data['f2'] == greedy) & (data['f3'] == doubleEdge)).nonzero()[0]

steps = np.unique(data['f6'][rows])
cpu_time = np.empty(np.size(steps),dtype = float)
q_min = np.empty(np.size(steps),dtype = float)
q_mean = np.empty(np.size(steps),dtype = float)
q_max = np.empty(np.size(steps),dtype = float)

for i in range(len(steps)):
    idx = (data['f6'] == steps[i]).nonzero()[0]
    cpu_time[i] = np.mean(data['f7'][idx])/1000
    q_min[i] = data['f8'][idx[0]]
    q_mean[i] = data['f9'][idx[0]]
    q_max[i]= data['f10'][idx[0]]
    
    
plt.figure()
plt.plot(steps,cpu_time)
plt.xlabel('Number of steps [-]')
plt.ylabel('CPU time [s]')
plt.xticks(range(0,np.max(steps)))


#%%
#from mpl_toolkits.mplot3d.art3d import Poly3DCollection
#fileName = '/stator_per_ffd_mod.su2'
##fileName = '/5x5x5_per_def.su2'
#[f_init,v_init,_,_,_,_,_] = getPlotData(fileName)  
#elem = v_init[f_init][76298]
#
#print(elem)
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1, projection='3d')
#
#verts = [np.array([elem[0],elem[1],elem[5],elem[4]]),
#      np.array([elem[3],elem[2],elem[6],elem[7]]),
#      np.array([elem[3],elem[0],elem[4],elem[7]]),
#      np.array([elem[2],elem[1],elem[5],elem[6]]),
#      np.array([elem[0],elem[1],elem[2],elem[3]]),
#      np.array([elem[4],elem[5],elem[6],elem[7]])]
#        
#pc = Poly3DCollection(verts, facecolors='red', edgecolor="black")
#ax.add_collection(pc)
#elem = v[f][76298]
#print(elem)
#verts = [np.array([elem[0],elem[1],elem[5],elem[4]]),
#      np.array([elem[3],elem[2],elem[6],elem[7]]),
#      np.array([elem[3],elem[0],elem[4],elem[7]]),
#      np.array([elem[2],elem[1],elem[5],elem[6]]),
#      np.array([elem[0],elem[1],elem[2],elem[3]]),
#      np.array([elem[4],elem[5],elem[6],elem[7]])]
#        
#pc = Poly3DCollection(verts, facecolors='blue', edgecolor="black")
#    
#ax.add_collection(pc)
#ax.set_xlim(0.290,0.303)
#ax.set_ylim(-0.065,-0.042)
#ax.set_zlim(0.0398,0.0408)
#plt.show()
#%%
if(0):
    #import numpy as np
    [f_init,v_init,elemType,bdryPnts_init,markerTags, nElemsMarks, FFD_pnts] = getPlotData("/stator_per_ffd_def.su2")
#    [f_init,v_init,elemType,bdryPnts_init,markerTags, nElemsMarks, FFD_pnts] = getPlotData("/OneraM6_fixed.su2")
#    plotTag= ["LOWER", "UPPER","LEFT","RIGHT","FRONT","BACK", "BLOCK"]
#    plotTag = ["BLADE", "HUB", "SHROUD"]
    plotTag = ["BLADE","HUB","INFLOW","OUTFLOW","PER1","PER2","SHROUD"]
#    plotTag = ["LOWER_SIDE" ,"UPPER_SIDE","TIP", "XNORMAL_FACES", "YNORMAL_FACES", "ZNORMAL_FACES", "SYMMETRY_FACE"]
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
if(0):
    import numpy as np
    import math
    import matplotlib.collections
    from colMap import colMap
    
    cmapMatlab = colMap()
#    fileNames = ['/OneraM6_fixed.su2']
    fileNames = ['/stator_per_ffd_def.su2']
#    fileNames = ['/5x5x5_per.su2']
    
#    [f_init,v_init,elemType,bdryPnts_init,_,_,_] = getPlotData('/OneraM6_fixed.su2')
    
    
    graphNames = ["",""]
    cutAxis = 0
    cutPlaneLoc = 0.27
    
    dims = [0,1,2]
    dims.remove(cutAxis)
    [f,v,elemType,_, markerTags, nElemsMarks,FFD_pnts] = getPlotData(fileNames[0]) 
       
    
    idx_lower = np.where(v[f][:,:,cutAxis] < cutPlaneLoc)
    idx_upper = np.where(v[f][:,:,cutAxis] >= cutPlaneLoc)
    idxCutElems = np.intersect1d(idx_lower[0],idx_upper[0])
    
    meshQuality = getMeshQual(fileNames[0])
    meshQual_cut = meshQuality[idxCutElems]
    

    print("Min mesh quality: \t", round(np.min(meshQual_cut),5)) 
    print("Max mesh quality: \t", round(np.max(meshQual_cut),5)) 
    print("Mean mesh quality: \t", round(np.mean(meshQual_cut),5))
    
    # quadriliteral
    verticesInd = np.array([[0,1],[1,2],[2,3],[3,0],[0,4],[1,5],[2,6],[3,7],[4,5],[5,6],[6,7],[7,4]])
    
    # tetrahidral
#    verticesInd = np.array([[0,1],[1,2],[2,0],[0,3],[1,3],[2,3]])
    
    v_cut = np.empty((len(idxCutElems),8,2),dtype = float)
    v_cut[:] = np.nan    
    
    for i in range(len(idxCutElems)):
        if(i%100 == 0):
            print(i)
#    for i in range(1):
        
        edges = v[f][idxCutElems[i]][verticesInd]   
        
#        cutLoc = np.empty((12,2), dtype = float)
#        cnt = 0
#        for ii in range(len(edges)):
#            if np.any(edges[ii][:,cutAxis] < cutPlaneLoc) and np.any(edges[ii][:,cutAxis] >= cutPlaneLoc):
#                delta = edges[ii][1,:] - edges[ii][0,:]
#                cutLoc[cnt,:] = edges[ii][0,dims] + delta[dims]*((cutPlaneLoc-edges[ii][0,cutAxis])/delta[cutAxis])
#                cnt = cnt+1
#                
        idx_1 = np.where(edges[:,:,cutAxis] < cutPlaneLoc)
        idx_2 = np.where(edges[:,:,cutAxis] >= cutPlaneLoc)
        idx_crossing = np.intersect1d(idx_1[0],idx_2[0])
        
        cutLoc = edges[idx_crossing][:,0,dims] + (edges[idx_crossing][:,1,dims]-edges[idx_crossing][:,0,dims])*((cutPlaneLoc-edges[idx_crossing][:,0,cutAxis])/(edges[idx_crossing][:,1,cutAxis]-edges[idx_crossing][:,0,cutAxis]))[:, np.newaxis]
        
        
        midpoint = np.mean(cutLoc,axis=0)
        
        
        
#        angles = np.empty(cnt,dtype=float)
#        for ii in range(cnt):
#            angles[ii] = math.atan2(cutLoc[ii,1]-midpoint[1],cutLoc[ii,0]-midpoint[0])
#    
#        print(angles)
        angles = np.arctan2(cutLoc[:,1]-midpoint[1],cutLoc[:,0]-midpoint[0])
        cutLoc = cutLoc[angles.argsort()]
        
        v_cut[i,0:np.shape(cutLoc)[0],:] = cutLoc
    #    print(cutLoc[0:cnt,:])
    
    #    plt.scatter(cutLoc[0:cnt,0],cutLoc[0:cnt,1], color= "blue") 
                
    #    for ii in range(len(v1)):
                
    #        print(v[f][idxCutElems[i]][v1,cutAxis], v[f][idxCutElems[i]][v2,cutAxis])
     

    # plotting stuff   
    colors = cmapMatlab(plt.Normalize(0,1)(meshQual_cut))
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

