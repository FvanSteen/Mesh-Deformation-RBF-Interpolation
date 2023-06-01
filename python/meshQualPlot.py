
import matplotlib.pyplot as plt
import matplotlib.collections
import numpy as np
from funs import getMeshQual, getPlotData
from colMap import colMap
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from sys import exit


plt.rcParams["font.family"] = "Arial"    
plt.rcParams["font.size"] = 12

def meshQualPlot2D(fileName, graphName, f, v, elemType):
    # importing the default matlab colormap from colMap.py
    plt.rcParams["figure.figsize"] = [6.4*2, 2*4.8]
    cmapMatlab = colMap()
    
    meshQual = getMeshQual(fileName)
       
    quad_idx = np.where(elemType == 9)
    tri_idx = np.where(elemType == 5)
    
    colors = cmapMatlab(plt.Normalize(0,1)(meshQual))

    pc_tri = matplotlib.collections.PolyCollection(v[f[tri_idx][:,0:3]],cmap=cmapMatlab, facecolors=colors, edgecolor="black",linewidth=0.05)
    pc_quad = matplotlib.collections.PolyCollection(v[f[quad_idx][:,0:4]],cmap=cmapMatlab,  facecolors=colors, edgecolor="black",linewidth=0.05)

    fig = plt.figure()              
    ax = fig.add_subplot(1,1,1)
    
    polys = ax.add_collection(pc_tri)
    polys = ax.add_collection(pc_quad)
    
#        pc.set_array(None)
    ax.autoscale()
    ax.set_aspect('equal')
    polys.set_clim(0,1)
    plt.colorbar(polys, ax=ax, shrink = .8)
    

    ax.title.set_text(graphName)
    plt.tight_layout()
    plt.show()

#    x = np.argwhere(np.isnan(meshQual))
#    print(v[f[x]])
    
    
def meshQualPlot3D(fileName, graphName, cutAxis, cutPlaneLoc, f, v, elemType):
    
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
    
    fig = plt.figure()  
    pc = matplotlib.collections.PolyCollection(v_cut,cmap=cmapMatlab, facecolors=colors, edgecolor="black",linewidth=0.1)
    ax = fig.add_subplot(1,1,1)
    polys = ax.add_collection(pc)
    
    v_cut2 = v_cut
    v_cut2[:,:,0] = v_cut2[:,:,0]-0.174533
    pc2 = matplotlib.collections.PolyCollection(v_cut2,cmap=cmapMatlab, facecolors=colors, edgecolor="black",linewidth=0.1)
    polys = ax.add_collection(pc2)
    
    
#    v_cut2[:,:,0] = v_cut2[:,:,0]-0.174533
#    pc2 = matplotlib.collections.PolyCollection(v_cut2,cmap=cmapMatlab, facecolors=colors, edgecolor="black",linewidth=0.1)
#    polys = ax.add_collection(pc2)
#    
#    v_cut2[:,:,0] = v_cut2[:,:,0]+3*0.174533
#    pc2 = matplotlib.collections.PolyCollection(v_cut2,cmap=cmapMatlab, facecolors=colors, edgecolor="black",linewidth=0.1)
#    polys = ax.add_collection(pc2)
#    
#    v_cut2[:,:,0] = v_cut2[:,:,0]+0.174533
#    pc2 = matplotlib.collections.PolyCollection(v_cut2,cmap=cmapMatlab, facecolors=colors, edgecolor="black",linewidth=0.1)
#    polys = ax.add_collection(pc2)
    
#    pc.set_array(None)
    ax.autoscale()  
#    ax.set_aspect('equal')
    plt.xlim([-0.25,-0.09])
    plt.ylim([0,0.06])
    polys.set_clim(0,1)
    plt.colorbar(polys, ax=ax,shrink=1)
    plt.xlabel('$\\theta$')
    plt.ylabel('z')
    plt.tight_layout()
    plt.show()

    
def nanElems(fileName,v,f):
    
    meshQuality = getMeshQual(fileName)
    idx = np.argwhere(np.isnan(meshQuality))
    idx_max = np.argwhere(meshQuality > 1)
    
    idx_min = np.argwhere(meshQuality < 0)
    print("idxmax", idx_max)
    print("idxmin", idx_min)
    print("idxnan", idx)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    for i in range(len(idx)):
        if(i%100 == 0):
            print(str(i) + "/" + str(len(idx)))
        elem = v[f][idx[i]][0]

        verts = [np.array([elem[0],elem[1],elem[5],elem[4]]),
              np.array([elem[3],elem[2],elem[6],elem[7]]),
              np.array([elem[3],elem[0],elem[4],elem[7]]),
              np.array([elem[2],elem[1],elem[5],elem[6]]),
              np.array([elem[0],elem[1],elem[2],elem[3]]),
              np.array([elem[4],elem[5],elem[6],elem[7]])]
        
        
        pc = Poly3DCollection(verts, facecolors='red', edgecolor="black",linewidth=0.5)
    
        ax.add_collection(pc)
        
    for i in range(len(idx_max)):
        if(i%100 == 0):
            print(str(i) + "/" + str(len(idx_max)))
        elem = v[f][idx_max[i]][0]
    
        verts = [np.array([elem[0],elem[1],elem[5],elem[4]]),
              np.array([elem[3],elem[2],elem[6],elem[7]]),
              np.array([elem[3],elem[0],elem[4],elem[7]]),
              np.array([elem[2],elem[1],elem[5],elem[6]]),
              np.array([elem[0],elem[1],elem[2],elem[3]]),
              np.array([elem[4],elem[5],elem[6],elem[7]])]
        
        
        pc = Poly3DCollection(verts, facecolors='yellow', edgecolor="black",linewidth=0.5)
    
        ax.add_collection(pc)
        
    for i in range(len(idx_min)):
        if(i%100 == 0):
            print(str(i) + "/" + str(len(idx_min)))
        elem = v[f][idx_min[i]][0]

        verts = [np.array([elem[0],elem[1],elem[5],elem[4]]),
              np.array([elem[3],elem[2],elem[6],elem[7]]),
              np.array([elem[3],elem[0],elem[4],elem[7]]),
              np.array([elem[2],elem[1],elem[5],elem[6]]),
              np.array([elem[0],elem[1],elem[2],elem[3]]),
              np.array([elem[4],elem[5],elem[6],elem[7]])]
        
        
        pc = Poly3DCollection(verts, facecolors='blue', edgecolor="black",linewidth=0.5)
    
        ax.add_collection(pc)
    
    ub = np.max(v,axis=0)
    lb = np.min(v,axis=0)
    ax.set_xlim(lb[0], ub[0])
    ax.set_ylim(lb[1], ub[1])
    ax.set_zlim(lb[2], ub[2])
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()
    

    
    