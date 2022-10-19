import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections
from colMap import colMap

def getMeshQuals(faces,vertices, alphas_0):
    
    [alphas,lambda_11, lambda_22] = getMeshQualParams(faces,vertices)
    
    tau = np.sum(alphas,axis=1)/np.sum(alphas_0, axis=1)
    f_size = np.amin(np.array([tau,1/tau]),axis=0)
    f_skew = 4/np.sum(np.sqrt(lambda_11*lambda_22)/alphas, axis=1)
    f_ss = np.sqrt(f_size)*f_skew
    return f_ss
    
def getMeshQualParams(faces, vertices):
    alphas = np.zeros(faces.shape)
    lambda_11 = np.zeros(faces.shape)
    lambda_22 = np.zeros(faces.shape)
    for i in range(0,faces.shape[0]):    
        elem = v[f[i]]
        for ii in range(0,faces.shape[1]):
            iip1 = (ii+1)%4
            iip3= (ii+3)%4
            A = np.array([[elem[iip1,0]-elem[ii,0], elem[iip3,0]-elem[ii,0]], [elem[iip1,1]-elem[ii,1], elem[iip3,1]-elem[ii,1]]])
            alphas[i,ii] = np.linalg.det(A)
            
            lambda_11[i,ii] = A[0,0]**2 + A[1,0]**2
            lambda_22[i,ii] = A[1,1]**2 + A[0,1]**2
            
    return [alphas,lambda_11,lambda_22]

def getPlotData(fileName, intBdryTag):
    fileObj = open(x+fileName, "r") 

    lines = fileObj.read().splitlines()
    fileObj.close()

    idx = 0
    for line in lines:
        if line.strip().startswith('NDIME= '):
            nDim = int(line[7:])
            nDimIdx = idx
        elif line.strip().startswith('NELEM= '):
            nElem = int(line[7:])
            nElemIdx = idx
        elif line.strip().startswith('NPOIN= '):
            nPnt = int(line[7:])
            nPntIdx = idx
        elif line.strip().startswith('MARKER_TAG= '+intBdryTag):
            intBdryIdx = idx+1
        
        idx += 1
    

    print("Dimensions:\t", nDim, "\nElements:\t",nElem, "\nPoints:\t\t", nPnt)
    print(nDimIdx,nElemIdx,nPntIdx)
    
    nLine = lines[intBdryIdx]
    nInElems = int(nLine[14:])
    f_in = np.empty([1,2*nInElems],dtype=int)
    
    for i in range(nInElems):
        lineData = lines[intBdryIdx+1+i].strip().split('\t')    
        f_in[0,2*i:2*i+2] = lineData[1:3]
    
    # Change the 3 and 4 here to be adjustable to the type of elements used in the mesh
    f = np.empty((nElem,4),dtype=int)
    for i in range(nElem):
        lineData = lines[i+nElemIdx+1].strip().split('\t')
        f[i,:] = lineData[1:5]
        
        
    
    v = np.empty((nPnt,nDim))
    for i in range(nPnt):
        lineData = lines[i+nPntIdx+1].strip().split('\t')
        v[i,:] = lineData[:nDim]
            
    return [f,v,f_in]

# importing the default matlab colormap from colMap.py
cmapMatlab = colMap();
 
# Setting directory to find .su2 files
os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\Classical RBF\\Meshes')
x = os.getcwd()

# Provide filenames of the meshes
fileNames = ['/5x5.su2', '/5x5_def.su2']
#fileNames = ['/TestMesh.su2', '/TestMesh_def.su2']
intBdryTag = "block"
#plt.close('all')
fig = plt.figure()


for i in range(0,len(fileNames)):
    ax = fig.add_subplot(1,len(fileNames),i+1)
    [f,v,f_in] = getPlotData(fileNames[i],intBdryTag)  
    
    if i==0:
        [alphas_0,_,_] = getMeshQualParams(f,v)
        
    meshQual = getMeshQuals(f,v,alphas_0)
    colors = cmapMatlab(plt.Normalize(0,1)(meshQual))
    pc = matplotlib.collections.PolyCollection(v[f],cmap=cmapMatlab, facecolors=colors, edgecolor="black",linewidth=0.5)
    pc_inner = matplotlib.collections.PolyCollection(v[f_in],cmap=cmapMatlab, facecolors="white", edgecolor="black",linewidth=0.5)
    polys = ax.add_collection(pc)
    polys = ax.add_collection(pc_inner)    
    
    pc.set_array(None)
    ax.autoscale()
    ax.set_aspect('equal')
    polys.set_clim(0,1)
    plt.colorbar(polys, ax=ax, shrink=0.44)
    ax.title.set_text('Initial Mesh')
    ax.set_ylim([-0.05, 1.05])


