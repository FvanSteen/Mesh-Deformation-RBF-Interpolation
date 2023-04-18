
import matplotlib.pyplot as plt
import matplotlib.collections
import numpy as np
from funs import getMeshQuals,getPlotData
from colMap import colMap


def meshQualPlot(fileNames,fNameInit, graphNames,alphas_0):
    # importing the default matlab colormap from colMap.py
    cmapMatlab = colMap()
    

    fig = plt.figure()  
    #fig.suptitle('NACA 0012',y=0.85)
        
    
    
    for i in range(0,len(fileNames)):
    #    ax = fig.add_subplot(1,len(fileNames),i+1)
        ax = fig.add_subplot(1,1,i+1)
        [f,v,elemType,_,_,_,_] = getPlotData(fileNames[i])  
#        print(v)
        
        coordsPolar = np.empty(np.shape(v))
        coordsPolar[:,0] = np.sqrt(v[:,0]**2+v[:,1]**2)
        coordsPolar[:,1] = np.arctan2(v[:,1],v[:,0])
        coordsPolar[:,1] += np.pi/6
        v2 = np.empty(np.shape(v))
        v2[:,0] = coordsPolar[:,0]*np.cos(coordsPolar[:,1])
        v2[:,1] = coordsPolar[:,0]*np.sin(coordsPolar[:,1])
        
        coordsPolar[:,1] += np.pi/6
        v3 = np.empty(np.shape(v))
        v3[:,0] = coordsPolar[:,0]*np.cos(coordsPolar[:,1])
        v3[:,1] = coordsPolar[:,0]*np.sin(coordsPolar[:,1])
        
        
        
#        coordsPolar[:,0] = np.sqrt(v[:,0]**2+v[:,1]**2)
#        print(coordsPolar)
        
        quadIdx = np.where(elemType == 9)
        if np.size(quadIdx) != 0:
            startQuadIdx = np.where(elemType == 9)[0][0]
        else:
            startQuadIdx = np.size(elemType)
        
        meshQual = getMeshQuals(f,v,alphas_0, elemType)
        
        print("Min mesh quality: \t", round(np.min(meshQual),5)) 
        print("Max mesh quality: \t", round(np.max(meshQual),5)) 
        print("Mean mesh quality: \t", round(np.mean(meshQual),5))
        
    #    colors1 = cmapMatlab(plt.Normalize(0,1)(meshQual[0:startQuadIdx]))
    #    
    #    colors2 = cmapMatlab(plt.Normalize(0,1)(meshQual[startQuadIdx:]))
        colors1 = cmapMatlab(plt.Normalize(0,1)(meshQual[0:startQuadIdx]))
        colors2 = cmapMatlab(plt.Normalize(0,1)(meshQual[startQuadIdx:]))
        
    #    colormap = plt.cm.coolwarm #or any other colormap
    #    cmap_reversed = matplotlib.cm.get_cmap('bwr_r')
    #    colors2 = cmap_reversed(plt.Normalize(0,1)(meshQual[startQuadIdx:]))
        
    #    v[:,0] = -v[:,0]
        pc = matplotlib.collections.PolyCollection(v[f[0:startQuadIdx,0:3]],cmap='seismic', facecolors=colors1, edgecolor="black",linewidth=0.1)
        pc2 = matplotlib.collections.PolyCollection(v[f[startQuadIdx:,0:4]],cmap=cmapMatlab,  facecolors=colors2, edgecolor="black",linewidth=0.1)
#        v[:,1] += 0.045
#        pc3 = matplotlib.collections.PolyCollection(v2[f[0:startQuadIdx,0:3]],cmap='seismic', facecolors=colors1, edgecolor="black",linewidth=0.1)
        pc4 = matplotlib.collections.PolyCollection(v2[f[startQuadIdx:,0:4]],cmap=cmapMatlab,  facecolors=colors2, edgecolor="black",linewidth=0.1)
        pc5 = matplotlib.collections.PolyCollection(v3[f[startQuadIdx:,0:4]],cmap=cmapMatlab,  facecolors=colors2, edgecolor="black",linewidth=0.1)
        
    #    pc3 = matplotlib.collections.PolyCollection(v[f[4365:4366]],cmap=cmapMatlab, facecolors='red', edgecolor="black",linewidth=0.25)
        
        polys = ax.add_collection(pc)
        polys = ax.add_collection(pc2)
#        polys = ax.add_collection(pc3)
#        polys = ax.add_collection(pc4)
#        polys = ax.add_collection(pc5)
        
    #    polys = ax.add_collection(pc3)
#        if(fileNames[i][0:10] == '/mesh_NACA'): 
#            ax.scatter(v[200][0],v[200][1],color='red') 
        pc.set_array(None)
        ax.autoscale()
        ax.set_aspect('equal')
        polys.set_clim(0,1)
        plt.colorbar(polys, ax=ax, shrink=1.0/(np.size(graphNames)-1))
        ax.title.set_text(graphNames[i])
#        ax.set_ylim([-0.031, -0.01])
#        ax.set_xlim([-0.140, -0.126])
    plt.show()
    x = np.argwhere(np.isnan(meshQual))
    print(v[f[x]])