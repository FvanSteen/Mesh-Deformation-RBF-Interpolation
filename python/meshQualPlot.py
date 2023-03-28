
import matplotlib.pyplot as plt
import matplotlib.collections
import numpy as np
from funs import getMeshQuals,getPlotData
from colMap import colMap


def meshQualPlot(meshPath, fileNames,fNameInit, alphas_0):
    # importing the default matlab colormap from colMap.py
    cmapMatlab = colMap()
    
    fig = plt.figure()   
    
    for i in range(0,len(fileNames)):
        ax = fig.add_subplot(1,1,i+1)
        [f,v,elemType,_,_,_,_] = getPlotData(meshPath + fileNames[i])  
        
        quadIdx = np.where(elemType == 9)
        if np.size(quadIdx) != 0:
            startQuadIdx = np.where(elemType == 9)[0][0]
        else:
            startQuadIdx = np.size(elemType)
        
        meshQual = getMeshQuals(f,v,alphas_0, elemType)
        
        print("Min mesh quality: \t", round(np.min(meshQual),5)) 
        print("Max mesh quality: \t", round(np.max(meshQual),5)) 
        print("Mean mesh quality: \t", round(np.mean(meshQual),5))
        

        colors1 = cmapMatlab(plt.Normalize(0,1)(meshQual[0:startQuadIdx]))
        colors2 = cmapMatlab(plt.Normalize(0,1)(meshQual[startQuadIdx:]))
        
        pc = matplotlib.collections.PolyCollection(v[f[0:startQuadIdx,0:3]],cmap='seismic', facecolors=colors1, edgecolor="black",linewidth=0.1)
        pc2 = matplotlib.collections.PolyCollection(v[f[startQuadIdx:,0:4]],cmap=cmapMatlab,  facecolors=colors2, edgecolor="black",linewidth=0.1)

        
        polys = ax.add_collection(pc)
        polys = ax.add_collection(pc2)

        pc.set_array(None)
        ax.autoscale()
        ax.set_aspect('equal')
        polys.set_clim(0,1)
        plt.colorbar(polys, ax=ax)

    plt.show()
#    x = np.argwhere(np.isnan(meshQual))
#    print(v[f[x]])