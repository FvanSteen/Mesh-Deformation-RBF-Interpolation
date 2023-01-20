import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections
from colMap import colMap
import math
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from getCrossSection import getCrossSection



def getVertices(elem):

    vertices = [np.array([elem[0],elem[1],elem[5],elem[4]]),
          np.array([elem[3],elem[2],elem[6],elem[7]]),
          np.array([elem[3],elem[0],elem[4],elem[7]]),
          np.array([elem[2],elem[1],elem[5],elem[6]]),
          np.array([elem[0],elem[1],elem[2],elem[3]]),
          np.array([elem[4],elem[5],elem[6],elem[7]])]
    return vertices 

def getMeshQuals(faces,vertices, alphas_0):
    
    [alphas,lambda_11, lambda_22,lambda_33] = getMeshQualParams(faces,vertices)
    
    tau = np.sum(alphas,axis=1)/np.sum(alphas_0, axis=1)
    f_size = np.amin(np.array([tau,1/tau]),axis=0)
    
#    f_skew = 4/np.sum(np.sqrt(lambda_11*lambda_22)/alphas, axis=1)
    
    f_skew = 8/np.sum(np.sqrt(lambda_11*lambda_22*lambda_33)/alphas, axis=1)
    f_ss = np.sqrt(f_size)*f_skew
    return f_ss
    
def getMeshQualParams(faces, vertices):
    alphas = np.zeros(faces.shape)
    lambda_11 = np.zeros(faces.shape)
    lambda_22 = np.zeros(faces.shape)
    lambda_33 = np.zeros(faces.shape)

    kp1 = [1,2,3,0,7,4,5,6]
    kp3 = [3,0,1,2,5,6,7,4]
    kp4 = [4,5,6,7,0,1,2,3]
    for i in range(0,faces.shape[0]): 
        elem = v[f[i]]

        for k in range(0,faces.shape[1]):

            A = np.array([[elem[kp1[k],0]-elem[k,0], elem[kp3[k],0]-elem[k,0], elem[kp4[k],0]-elem[k,0]],
                          [elem[kp1[k],1]-elem[k,1], elem[kp3[k],1]-elem[k,1], elem[kp4[k],1]-elem[k,1]],
                          [elem[kp1[k],2]-elem[k,2], elem[kp3[k],2]-elem[k,2], elem[kp4[k],2]-elem[k,2]]])

            alphas[i,k] = np.linalg.det(A)

            tensor = np.matmul(np.transpose(A),A)
            
            lambda_11[i,k] = tensor[0,0]
            lambda_22[i,k] = tensor[1,1]
            lambda_33[i,k] = tensor[2,2]
    return [alphas,lambda_11,lambda_22,lambda_33]


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
        elif line.strip().startswith('NMARK='):
            bdryIdx = idx+3
        idx += 1
        
    

    print("Dimensions:\t", nDim, "\nElements:\t",nElem, "\nPoints:\t\t", nPnt)
    print(nDimIdx,nElemIdx,nPntIdx)
#    
    nLine = lines[intBdryIdx]
    nInElems = int(nLine[14:])
    f_in = []   
    
    for i in range(nInElems):
        lineData = lines[intBdryIdx+1+i].strip().split('\t')    
        
#        f_in[0,2*i:2*i+2] = lineData[1:3]
        for j in lineData[1:5]:
            f_in.append(int(j))
        
    # Change the 3 and 4 here to be adjustable to the type of elements used in the mesh
    
    if elementType == 5:
        nNodesElem = 3
    elif elementType == 9:
        nNodesElem = 4
    elif elementType == 12:
        nNodesElem = 8
#    nElem = 1
    f = np.empty((nElem,nNodesElem),dtype=int)
    for i in range(nElem):
        lineData = lines[i+nElemIdx+1].strip().split('\t')
        f[i,:] = lineData[1:nNodesElem+1]
        
        
    
    v = np.empty((nPnt,nDim))
    for i in range(nPnt):
        lineData = lines[i+nPntIdx+1].strip().split('\t')
        v[i,:] = lineData[:nDim]
            
    f_in = np.unique(f_in)
    
    bdryNodes = []
    for i in range(bdryIdx,intBdryIdx):
        lineData = lines[i].strip().split('\t')   
        if lineData[0].isdigit():
            for j in lineData[1:-1]:
                bdryNodes.append(int(j))
                
    bdryNodes = np.unique(bdryNodes)
#    print(bdryNodes)
    
    
    
    return [f,v,f_in,bdryNodes]

# importing the default matlab colormap from colMap.py
cmapMatlab = colMap();
 
# Setting directory to find .su2 files
os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\Classical RBF\\Meshes')
#os.chdir('c:\\Users\\floyd\\eclipse-workspace\\Thesis\\createMeshFiles')



# Provide filenames of the meshes
#fileNames = ['/25x25mesh.su2', '/25x25mesh_def.su2']
#fileNames = ['/TestMesh.su2', '/TestMesh_def.su2']
intBdryTag = "BLOCK"
#fileNames = ["/25x25x5_def_rbf.su2"] 
fileNames = ["/5x5x5.su2"]
gt = 'RBF'

# Obtaining mesh quality parameters from the initial mesh
[f,v,f_in,_] = getPlotData("/5x5x5.su2",intBdryTag)
[alphas_0,_,_,_] = getMeshQualParams(f,v)
# Element identifiers according to the VTK format:
# Triangle:         5
# Quadrilateral:    9
# Hexahedral:       12
# 
elementType = 12 
[f,v,f_in,bdryNodes] = getPlotData(fileNames[0],elementType,intBdryTag)
meshQual = getMeshQuals(f,v,alphas_0)


z_cut = 0.5
dim = 2
[data, ver_in,mq] = getCrossSection(f,v,z_cut,dim,meshQual,f_in)

colors = cmapMatlab(plt.Normalize(0,1)(mq))

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
pc = matplotlib.collections.PolyCollection(data,cmap=cmapMatlab, facecolors=colors, edgecolor='black')
pc_in = matplotlib.collections.PolyCollection(ver_in, facecolors='white', edgecolor='black')
polys = ax.add_collection(pc)
polys = ax.add_collection(pc_in)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_xlim(0,1)
ax.set_ylim(0,1)
ax.set_aspect('equal')
ax.title.set_text(gt + '\nCross section at z='+str(z_cut))
pc.set_array(None)
polys.set_clim(0,1)
plt.colorbar(polys, ax=ax, shrink=0.44)
plt.show()
#%%


colors = cmapMatlab(plt.Normalize(0,1)(meshQual))

#fig = plt.figure()
#ax = fig.add_subplot(1,1,1,projection='3d')
##for i in range(len(f)):
#idx = list(range(2,len(f),5))
#for i in idx:
#    verts = [np.array([v[f][i][0],v[f][i][1],v[f][i][5],v[f][i][4]]),
#          np.array([v[f][i][3],v[f][i][2],v[f][i][6],v[f][i][7]]),
#          np.array([v[f][i][3],v[f][i][0],v[f][i][4],v[f][i][7]]),
#          np.array([v[f][i][2],v[f][i][1],v[f][i][5],v[f][i][6]]),
#          np.array([v[f][i][0],v[f][i][1],v[f][i][2],v[f][i][3]]),
#          np.array([v[f][i][4],v[f][i][5],v[f][i][6],v[f][i][7]])]
#    pc = Poly3DCollection(verts, cmap=cmapMatlab, facecolors=colors[i], edgecolor="black",linewidth=0.5)
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
#ax.azim = -125
#ax.elev = 23
#ax.title.set_text(gt)
#%% plotting of internal boundary elements
idxIntBdryElems = []
for i in range(len(f)):
    if set(f[i]).issubset(f_in):
        idxIntBdryElems.append(i)
   

        

elem = v[f][idxIntBdryElems][0]
verts = getVertices(elem)

fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection='3d')
pc = Poly3DCollection(verts, cmap=cmapMatlab, facecolors='blue', edgecolor="black",linewidth=0.5)
polys = ax.add_collection(pc)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

ax.set_xlim(0,1)
ax.set_ylim(0,1)
ax.set_zlim(0,1)
ax.autoscale()

#%% plotting of outer boundaries

idxBdryElems = []
for i in range(len(f)):
    if np.any(np.in1d(f[i], bdryNodes)):
        idxBdryElems.append(i)
    

fig = plt.figure()        
ax = fig.add_subplot(1,1,1, projection='3d')
for i in range(len(idxBdryElems)):
    elem = v[f][idxBdryElems][i]
    verts = getVertices(elem)
    pc = Poly3DCollection(verts, cmap=cmapMatlab, facecolors=colors[idxBdryElems[i]], edgecolor="black",linewidth=0.5)
    polys = ax.add_collection(pc)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

ax.set_xlim(0,1)
ax.set_ylim(0,1)
ax.set_zlim(0,1)
ax.autoscale()
ax.title.set_text(gt+ '\n3D view')
polys.set_clim(0,1)
plt.colorbar(polys, ax=ax, shrink=0.44)
ax.azim = -125
ax.elev = 23
