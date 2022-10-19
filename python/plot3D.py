import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections
from colMap import colMap
import math
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from getCrossSection import getCrossSection
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
    x = os.getcwd()
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
#    
    nLine = lines[intBdryIdx]
    nInElems = int(nLine[14:])
    f_in = np.empty([1,2*nInElems],dtype=int)
#    
#    for i in range(nInElems):
#        lineData = lines[intBdryIdx+1+i].strip().split('\t')    
#        f_in[0,2*i:2*i+2] = lineData[1:3]
    
    # Change the 3 and 4 here to be adjustable to the type of elements used in the mesh
#    nElem = 1
    f = np.empty((nElem,8),dtype=int)
    for i in range(nElem):
        lineData = lines[i+nElemIdx+1].strip().split('\t')
        f[i,:] = lineData[1:9]
        
        
    
    v = np.empty((nPnt,nDim))
    for i in range(nPnt):
        lineData = lines[i+nPntIdx+1].strip().split('\t')
        v[i,:] = lineData[:nDim]
            
        
    return [f,v,f_in]

# importing the default matlab colormap from colMap.py
cmapMatlab = colMap();
 
# Setting directory to find .su2 files
os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\Classical RBF\\Meshes')
#os.chdir('c:\\Users\\floyd\\eclipse-workspace\\Thesis\\createMeshFiles')



# Provide filenames of the meshes
#fileNames = ['/25x25mesh.su2', '/25x25mesh_def.su2']
#fileNames = ['/TestMesh.su2', '/TestMesh_def.su2']
intBdryTag = "BLOCK"
fileNames = ["/15x15x5_def.su2"] 

[f,v,f_in] = getPlotData(fileNames[0],intBdryTag)



#verts = [np.array([[0,0,0], [1,0,0], [1,1,1], [0,1,1]])]

z_cut = .5
data = getCrossSection(f,v,z_cut)


fig = plt.figure()
ax = fig.add_subplot(1,1,1)
pc = matplotlib.collections.PolyCollection(data, facecolor='b', edgecolor='black')
polys = ax.add_collection(pc)
ax.set_xlim(0,1)
ax.set_ylim(0,1)
ax.set_aspect('equal')
plt.show()

#
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1,projection='3d')
##for i in range(len(f)):
#for i in range(50,75):
#    verts = [np.array([v[f][i][0],v[f][i][1],v[f][i][5],v[f][i][4]]),
#          np.array([v[f][i][3],v[f][i][2],v[f][i][6],v[f][i][7]]),
#          np.array([v[f][i][3],v[f][i][0],v[f][i][4],v[f][i][7]]),
#          np.array([v[f][i][2],v[f][i][1],v[f][i][5],v[f][i][6]]),
#          np.array([v[f][i][0],v[f][i][1],v[f][i][2],v[f][i][3]]),
#          np.array([v[f][i][4],v[f][i][5],v[f][i][6],v[f][i][7]])]
#    pc = Poly3DCollection(verts, facecolors="red", edgecolor="black",linewidth=0.5)
#    polys = ax.add_collection3d(pc)
#
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.set_zlabel('z')
##
###ax.azim = 90
###ax.elev = 90
#ax.set_xlim(0,1)
#ax.set_ylim(0,1)
#ax.set_zlim(0,1)







