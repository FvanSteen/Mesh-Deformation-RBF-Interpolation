# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 11:34:07 2023

@author: floyd
"""


import matplotlib.pyplot as plt
import matplotlib.collections
import numpy as np
from funs import getMeshQual, getPlotData
from colMap import colMap
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from sys import exit
from coordTransform import coordTransform
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from funs import getMeshQualParams,getPlotData, getMeshQualParams3D, getMeshQuals3D, getMeshQual, getMeshQuals
import os
figPath = os.path.dirname(os.path.abspath(__file__)) + "/figs/"
os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\Meshes')

fileName = '/stator_per_ffd_elastic_deform.su2'
#fileName = '/25x25x25.su2'
[f,v,elemType,bdryPnts,markerTags, nElemsMarks, FFD_pnts] = getPlotData(fileName)  
v = coordTransform.toCylindrical(v)
graphName = ''
cutAxis = 0
cutPlaneLoc = 0.2512

plt.rcParams["font.size"] = 14

#    plt.savefig(figPath + "stator_r=0.2512_ps_double.png", dpi=300,bbox_inches='tight')

plt.rcParams["font.family"] = "Arial"    
    
cmapMatlab = colMap()
   
dims = [0,1,2]
dims.remove(cutAxis)
   
meshQuality = getMeshQual(fileName)

idx_lower = np.where(v[f][:,:,cutAxis] < cutPlaneLoc)
idx_upper = np.where(v[f][:,:,cutAxis] >= cutPlaneLoc)
idxCutElems = np.intersect1d(idx_lower[0],idx_upper[0])

meshQual_cut = meshQuality[idxCutElems]

print("\nMesh quality parameters of the cross-section:" )
print("Min mesh quality: \t", round(np.min(meshQual_cut),5)) 
print("Max mesh quality: \t", round(np.max(meshQual_cut),5)) 
print("Mean mesh quality: \t", round(np.mean(meshQual_cut),5), "\n")

# quadriliteral
verticesInd = np.array([[0,1],[1,2],[2,3],[3,0],[0,4],[1,5],[2,6],[3,7],[4,5],[5,6],[6,7],[7,4]])

# tetrahidral
#    verticesInd = np.array([[0,1],[1,2],[2,0],[0,3],[1,3],[2,3]])

v_cut = np.empty((len(idxCutElems),8,2),dtype = float)
v_cut[:] = np.nan    

for i in range(len(idxCutElems)):
    if(i%100 == 0):
        print(str(i) + "/" + str(len(idxCutElems)))
    
    edges = v[f][idxCutElems[i]][verticesInd]   
    
    idx_1 = np.where(edges[:,:,cutAxis] < cutPlaneLoc)
    idx_2 = np.where(edges[:,:,cutAxis] >= cutPlaneLoc)
    idx_crossing = np.intersect1d(idx_1[0],idx_2[0])
    
    cutLoc = edges[idx_crossing][:,0,dims] + (edges[idx_crossing][:,1,dims]-edges[idx_crossing][:,0,dims])*((cutPlaneLoc-edges[idx_crossing][:,0,cutAxis])/(edges[idx_crossing][:,1,cutAxis]-edges[idx_crossing][:,0,cutAxis]))[:, np.newaxis]
    
    midpoint = np.mean(cutLoc,axis=0)
    
    angles = np.arctan2(cutLoc[:,1]-midpoint[1],cutLoc[:,0]-midpoint[0])
    cutLoc = cutLoc[angles.argsort()]
    
    v_cut[i,0:np.shape(cutLoc)[0],:] = cutLoc
    
print(str(len(idxCutElems)) + "/" + str(len(idxCutElems)))

colors = cmapMatlab(plt.Normalize(0,1)(meshQual_cut))

fig, ax = plt.subplots(1, 2, sharey=True)
plt.subplots_adjust(left=None, right=None, wspace=0.033)
pc = matplotlib.collections.PolyCollection(v_cut,cmap=cmapMatlab, facecolors=colors, edgecolor="black",linewidth=0.1)
polys = ax[0].add_collection(pc)
ax[0].axis('off')
v_cut2 = v_cut
v_cut2[:,:,0] = v_cut2[:,:,0]-0.174533
pc3 = matplotlib.collections.PolyCollection(v_cut2,cmap=cmapMatlab, facecolors=colors, edgecolor="black",linewidth=0.1)
polys = ax[0].add_collection(pc3)
ax[0].set_xlim(-0.25,-0.09)
ax[0].set_ylim(0,0.06)
fileName = '/25x25x25_def.su2'
fileName = '/stator_per_ffd_ps.su2'
[f,v,elemType,bdryPnts,markerTags, nElemsMarks, FFD_pnts] = getPlotData(fileName) 
v = coordTransform.toCylindrical(v)
meshQuality = getMeshQual(fileName)

idx_lower = np.where(v[f][:,:,cutAxis] < cutPlaneLoc)
idx_upper = np.where(v[f][:,:,cutAxis] >= cutPlaneLoc)
idxCutElems = np.intersect1d(idx_lower[0],idx_upper[0])

meshQual_cut = meshQuality[idxCutElems]

print("\nMesh quality parameters of the cross-section:" )
print("Min mesh quality: \t", round(np.min(meshQual_cut),5)) 
print("Max mesh quality: \t", round(np.max(meshQual_cut),5)) 
print("Mean mesh quality: \t", round(np.mean(meshQual_cut),5), "\n")

# quadriliteral
verticesInd = np.array([[0,1],[1,2],[2,3],[3,0],[0,4],[1,5],[2,6],[3,7],[4,5],[5,6],[6,7],[7,4]])

# tetrahidral
#    verticesInd = np.array([[0,1],[1,2],[2,0],[0,3],[1,3],[2,3]])

v_cut = np.empty((len(idxCutElems),8,2),dtype = float)
v_cut[:] = np.nan    

for i in range(len(idxCutElems)):
    if(i%100 == 0):
        print(str(i) + "/" + str(len(idxCutElems)))
    
    edges = v[f][idxCutElems[i]][verticesInd]   
    
    idx_1 = np.where(edges[:,:,cutAxis] < cutPlaneLoc)
    idx_2 = np.where(edges[:,:,cutAxis] >= cutPlaneLoc)
    idx_crossing = np.intersect1d(idx_1[0],idx_2[0])
    
    cutLoc = edges[idx_crossing][:,0,dims] + (edges[idx_crossing][:,1,dims]-edges[idx_crossing][:,0,dims])*((cutPlaneLoc-edges[idx_crossing][:,0,cutAxis])/(edges[idx_crossing][:,1,cutAxis]-edges[idx_crossing][:,0,cutAxis]))[:, np.newaxis]
    
    midpoint = np.mean(cutLoc,axis=0)
    
    angles = np.arctan2(cutLoc[:,1]-midpoint[1],cutLoc[:,0]-midpoint[0])
    cutLoc = cutLoc[angles.argsort()]
    
    v_cut[i,0:np.shape(cutLoc)[0],:] = cutLoc
    
print(str(len(idxCutElems)) + "/" + str(len(idxCutElems)))

colors = cmapMatlab(plt.Normalize(0,1)(meshQual_cut))





pc2 = matplotlib.collections.PolyCollection(v_cut,cmap=cmapMatlab, facecolors=colors, edgecolor="black",linewidth=0.1)
polys2 = ax[1].add_collection(pc2)
ax[1].axis('off')
v_cut2 = v_cut
v_cut2[:,:,0] = v_cut2[:,:,0]-0.174533
pc4 = matplotlib.collections.PolyCollection(v_cut2,cmap=cmapMatlab, facecolors=colors, edgecolor="black",linewidth=0.1)
polys2 = ax[1].add_collection(pc4)


#v_cut2[:,:,0] = v_cut2[:,:,0]-0.174533
#pc2 = matplotlib.collections.PolyCollection(v_cut2,cmap=cmapMatlab, facecolors=colors, edgecolor="black",linewidth=0.1)
#polys = ax.add_collection(pc2)
#    
#    v_cut2[:,:,0] = v_cut2[:,:,0]+3*0.174533
#    pc2 = matplotlib.collections.PolyCollection(v_cut2,cmap=cmapMatlab, facecolors=colors, edgecolor="black",linewidth=0.1)
#    polys = ax.add_collection(pc2)
#    
#    v_cut2[:,:,0] = v_cut2[:,:,0]+0.174533
#    pc2 = matplotlib.collections.PolyCollection(v_cut2,cmap=cmapMatlab, facecolors=colors, edgecolor="black",linewidth=0.1)
#    polys = ax.add_collection(pc2)

#    pc.set_array(None)
#ax.autoscale()  
#    ax.set_aspect('equal')
ax[1].set_xlim(-0.25,-0.09)
ax[1].set_ylim(0,0.06)
#plt.ylim([0,0.06])
#    plt.xlim([-0.195,-0.170])
#    plt.ylim([0.0250,0.0450])



#plt.tight_layout()
plt.savefig(figPath + "fp.png", format = 'png', dpi=600,bbox_inches='tight', transparant=True)
plt.show()

    
