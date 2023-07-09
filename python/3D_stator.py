import os
from funs import getPlotData
from bdryScatter import bdryScatterFuns
os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\Meshes')



fileNames = ["\\stator_per_ffd.su2"]
[f_init,v_init,elemType,bdryPnts_init, markerTags, nElemsMarks,FFD_pnts] = getPlotData(fileNames[0])
#plotTag = ["BLADE","HUB","INFLOW","OUTFLOW","PER1","PER2","SHROUD"]
#plotTag = ["HUB", "PER1","PER2"]
#bdryScatterFuns.bdryScatter3D(v_init, bdryPnts_init, markerTags, nElemsMarks, plotTag, FFD_pnts)

if(0): 
#%%

    fileNames = ["\\stator_per_ffd_deform.su2"]
    [f_init,v_init,elemType,bdryPnts_init, markerTags, nElemsMarks,FFD_pnts] = getPlotData(fileNames[0])
    plotTag = ["BLADE","HUB","INFLOW","OUTFLOW","PER1","PER2","SHROUD"]
    bdryScatterFuns.bdryScatter3D(v_init, bdryPnts_init, markerTags, nElemsMarks, plotTag, FFD_pnts)
    
    
    #%% getting the displacement of the blade
    import numpy as np
    import matplotlib.pyplot as plt
    fileNames = ["\\stator_per_ffd_mod.su2"]
    [f_init,v_init,elemType,bdryPnts_init, markerTags, nElemsMarks,FFD_pnts] = getPlotData(fileNames[0])
    fileNames = ["\\stator_per_ffd_elastic_deform.su2"]
    [_,v,elemType,_, markerTags, nElemsMarks,FFD_pnts] = getPlotData(fileNames[0])
    
    blade_marker = "BLADE"
    
    idx = markerTags.index(blade_marker)                        
    bdryPntsPlot =  np.unique(bdryPnts_init[sum(nElemsMarks[0:idx]):sum(nElemsMarks[0:idx+1]),:])
#    bdryPntsPlot = np.unique(bdryPnts_init)
    
    
    delta = v[bdryPntsPlot] - v_init[bdryPntsPlot]
    #%%
    def getDefFile(defFileName, index, disp):
        x = os.getcwd()
        o = open(x + defFileName, "w")
        
        
        
        for idx,c1,c2,c3 in zip(np.transpose(index),disp[:,0],disp[:,1],disp[:,2]):
            o.write("{0}\t{1}\t{2}\t{3}".format(idx, c1,c2,c3) + "\n")
        o.close()
    
    
    
    os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\defs')
    defFileName = "\\stator_per_ffd_deformation_all.txt"
#    
    getDefFile(defFileName, bdryPntsPlot, delta)
    
    #%%
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    
    ax.scatter(v_init[bdryPntsPlot][:,0],v_init[bdryPntsPlot][:,1],v_init[bdryPntsPlot][:,2],marker = "1")
    ax.scatter(v[bdryPntsPlot][:,0],v[bdryPntsPlot][:,1],v[bdryPntsPlot][:,2],marker = "1")
    plt.show()
    
#%%
    
    
if(0):
    fileName = "\\stator_per_ffd.su2"
    fileObj = open('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\Meshes'+fileName, "r") 
    lines = fileObj.read().splitlines()
    fileObj.close()

    markIdx = []
    nElems = []
    idx = 0

    for line in lines:
        if line.strip().startswith('MARKER_TAG= '):
            if(line.strip().find("SEND_RECEIVE") != -1):
                markIdx.append(idx)
                
        idx += 1
        
        
    
    for index in markIdx:
        nElems.append(int(lines[index+1][14:].split()[0]))
    
    nodes = np.empty((max(nElems),len(markIdx)), dtype = int)
    for i in range(len(nElems)):
        for j in range(nElems[i]):
            nodes[j,i] = int(lines[markIdx[i]+3+j].strip().split()[1])
            
#    bdryPntsPlot = nodes[:,0]
#    bdryPntsPlot2 = nodes[int(nElems[1]/2):,1]    
#    fig = plt.figure()
#    ax = fig.add_subplot(projection='3d')
#    ax.scatter(v[bdryPntsPlot][:,0],v[bdryPntsPlot][:,1],v[bdryPntsPlot][:,2],marker = "1", alpha = .2)
#    ax.scatter(v[bdryPntsPlot2][:,0],v[bdryPntsPlot2][:,1],v[bdryPntsPlot2][:,2],marker = "1", alpha = .2)
    

#%%
#    to be removed indices:
    idx_gone = nodes[int(nElems[1]/2):,1]
    
    idx_f_gone = []

    for i in range(len(f_init)):
        common = np.isin(f_init[i], idx_gone)
        if np.any(common) == True:
            idx_f_gone.append(i)
    f_new = np.delete(f_init,idx_f_gone,axis=0)        
    
    #%%
    markerL = []
    for marker in range(len(markerTags)):
        start_idx = int(np.sum(nElemsMarks[:marker]))
        cnt = 0
        for elem in range(nElemsMarks[marker]):            
            if np.any(np.isin(idx_gone,bdryPnts_init[start_idx+elem])):
                cnt += 1
        print(cnt)
        markerL.append(nElemsMarks[marker]-cnt)
    
    print(markerL)
    
    
#%%           
    o = open('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\Meshes' + "\\stator_per_ffd_test.su2", "w")
    i = 0
    markerCnt = 0
    while(i < len(lines)):
        if lines[i].strip().startswith('NELEM='):
            o.write('NELEM= ' + str(len(f_init)-len(idx_f_gone)) + "\n")
            elem_start_idx = i+1
            for ii in range(len(f_init)):
                common = np.isin(idx_gone,f_init[ii])
                if np.all(common == False):
                    o.write(lines[ii+elem_start_idx] + '\n')
            i += len(f_init)+1
        elif lines[i].strip().startswith('NPOIN='):
            o.write('NPOIN= ' + str(len(v_init)-len(idx_gone)) + "\n")
            point_start_idx = i+1
            for ii in range(len(v_init)):
                if np.isin(ii, idx_gone) == False:
                    o.write(lines[ii + point_start_idx] + "\n")
            i += len(v_init)+1
            
        elif lines[i].strip().startswith('NMARK='): 
            o.write('NMARK= 7\n')
            i+=1
        elif lines[i].strip().startswith('MARKER_ELEMS='):
            o.write('MARKER_ELEMS= ' + str(markerL[markerCnt]) + "\n")
            bdry_start_idx = int(np.sum(nElemsMarks[:marker]))
            start_idx = i+1
            tel = 0
            print(nElemsMarks[markerCnt])
            for ii in range(nElemsMarks[markerCnt]):
                if np.all(np.isin(bdryPnts_init[bdry_start_idx+ii],idx_gone) == False):
                    o.write(lines[start_idx + ii] + "\n")
                else: 
                    tel += 1
                    print(bdryPnts_init[bdry_start_idx+ii])
            print("COUNT: ",tel)
#            if markerCnt == 3:
#                fdsafd
            i+=nElemsMarks[markerCnt]+1
            markerCnt += 1
            
        elif lines[i].strip().startswith('MARKER_TAG= SEND_RECEIVE'):
            n = int(lines[i+1][14:].split()[0])+3
            i += n
        else:
            o.write(lines[i] + "\n")
            i += 1
        
    o.close()
        
        
    

        
#            [f_init,v_init,elemType,bdryPnts_init,markerTags, nElemsMarks, FFD_pnts] = getPlotData('/stator_per_ffd_base.su2')
#        [f,v,elemType,bdryPnts,markerTags, nElemsMarks, FFD_pnts] = getPlotData('/stator_per_ffd_elastic_deform.su2')
#        blade_marker = "BLADE"
#        idx = markerTags.index(blade_marker)                        
#        bdryPntsPlot =  np.unique(bdryPnts_init[sum(nElemsMarks[0:idx]):sum(nElemsMarks[0:idx+1]),:])
#        ax.scatter(v_init[bdryPntsPlot][:,0],v_init[bdryPntsPlot][:,1],v_init[bdryPntsPlot][:,2], marker = "1", alpha = .4)
#        ax.scatter(v[bdryPntsPlot][:,0],v[bdryPntsPlot][:,1],v[bdryPntsPlot][:,2], marker = "1", alpha = .4)