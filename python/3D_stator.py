import os
from funs import getPlotData
from bdryScatter import bdryScatterFuns
os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\Meshes')



fileNames = ["\\stator_per_ffd.su2"]
[f_init,v_init,elemType,bdryPnts_init, markerTags, nElemsMarks,FFD_pnts] = getPlotData(fileNames[0])
plotTag = ["BLADE","HUB","INFLOW","OUTFLOW","PER1","PER2","SHROUD"]
#plotTag = ["HUB", "PER1","PER2"]
bdryScatterFuns.bdryScatter3D(v_init, bdryPnts_init, markerTags, nElemsMarks, plotTag, FFD_pnts)

if(0):
#%%
    
    fileNames = ["\\stator_per_ffd_deform.su2"]
    [f_init,v_init,elemType,bdryPnts_init, markerTags, nElemsMarks,FFD_pnts] = getPlotData(fileNames[0])
    plotTag = ["BLADE","HUB","INFLOW","OUTFLOW","PER1","PER2","SHROUD"]
    bdryScatterFuns.bdryScatter3D(v_init, bdryPnts_init, markerTags, nElemsMarks, plotTag, FFD_pnts)
    
    
    #%% getting the displacement of the blade
    import numpy as np
    import matplotlib.pyplot as plt
    fileNames = ["\\stator_per_ffd.su2"]
    [f_init,v_init,elemType,bdryPnts_init, markerTags, nElemsMarks,FFD_pnts] = getPlotData(fileNames[0])
    fileNames = ["\\stator_per_ffd_deform.su2"]
    [_,v,elemType,_, markerTags, nElemsMarks,FFD_pnts] = getPlotData(fileNames[0])
    
    blade_marker = "BLADE"
    
    idx = markerTags.index(blade_marker)        
                    
    bdryPntsPlot =  np.unique(bdryPnts_init[sum(nElemsMarks[0:idx]):sum(nElemsMarks[0:idx+1]),:])
    delta = v[bdryPntsPlot] - v_init[bdryPntsPlot]
    #%%
    def getDefFile(defFileName, index, disp):
        x = os.getcwd()
        o = open(x + defFileName, "w")
        
        
        
        for idx,c1,c2,c3 in zip(np.transpose(index),disp[:,0],disp[:,1],disp[:,2]):
            o.write("{0}\t{1}\t{2}\t{3}".format(idx, c1,c2,c3) + "\n")
        o.close()
    
    
    
    os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\defs')
    defFileName = "\\stator_per_ffd_deformation.txt"
    
    getDefFile(defFileName, bdryPntsPlot, delta)
    
    #%%
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    
    ax.scatter(v_init[bdryPntsPlot][:,0],v_init[bdryPntsPlot][:,1],v_init[bdryPntsPlot][:,2],marker = "1")
    ax.scatter(v[bdryPntsPlot][:,0],v[bdryPntsPlot][:,1],v[bdryPntsPlot][:,2],marker = "1")
    plt.show()
    
    
    
    
        