import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.collections
#def plot_grid(x,y, ax=None, **kwargs):
#    ax = ax or plt.gca()
#    segs1 = np.stack((x,y), axis=2)
#    segs2 = segs1.transpose(1,0,2)
#    ax.add_collection(LineCollection(segs1, **kwargs))
#    ax.add_collection(LineCollection(segs2, **kwargs))
#    ax.autoscale()


#os.chdir('c:\\Users\\floyd\\OneDrive\\Documenten\\Thesis\\Mesh-Deformation-RBF-Interpolation\\Classical RBF' )
os.chdir('../Classical RBF' )
x = os.getcwd()
fileName = '/mesh_NACA0012_inv.su2'
#fileName = '/TestMesh.su2'

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
    idx += 1
    

print("Dimensions:\t", nDim, "\nElements:\t",nElem, "\nPoints:\t\t", nPnt)
print(nDimIdx,nElemIdx,nPntIdx)


# Change the 3 and 4 here to be adjustable to the type of elements used in the mesh
f = np.empty((nElem,3),dtype=int)
for i in range(nElem):
    lineData = lines[i+nElemIdx+1].strip().split('\t')
    f[i,:] = lineData[1:4]
    
    
    
v = np.empty((nPnt,nDim))
for i in range(nPnt):
    lineData = lines[i+nPntIdx+1].strip().split('\t')
    v[i,:] = lineData[:nDim]
    
C = np.ones((1,nElem))


fig = plt.figure()
#ax = fig.add_subplot(111,projection="3d")
ax = plt.gca()
 

norm = plt.Normalize(C[0].min(), C[0].max())
colors = plt.cm.viridis(norm(C[0]))
pc = matplotlib.collections.PolyCollection(v[f],facecolor = 'white', edgecolor="black")
#ax.add_collection3d(pc)
ax.add_collection(pc)
ax.autoscale()

#%%
coords = v[f[1,:]]
#print(coords)
A1 = coords[1:3,:]
#print()
#print(A1)
#print()
A2 = coords[0:3:2,:]
#print(A2)

#print()
A3 = coords[0:2,:]
#print(A3)

a1 = np.linalg.det(A1)
a2 = np.linalg.det(A2)
a3 = np.linalg.det(A3)

tens1 = np.transpose(A1)*A1

tens2 = np.transpose(A2)*A2
tens3 = np.transpose(A3)*A3
#print (tens1,'\n','\n',tens2,'\n','\n',tens3)

print(np.sqrt(3)*a1/(tens1[0,0] + tens1[1,1] - tens1[0,1]))
print(np.sqrt(3)*a2/(tens2[0,0] + tens2[1,1] - tens2[0,1]))
print(np.sqrt(3)*a3/(tens3[0,0] + tens3[1,1] - tens3[0,1]))






#%%

nDimsIdx = [lines.index(l) for l in lines if l.startswith('NDIME= ')]
nDimsLine = lines[nDimsIdx[0]]
nDims = int(nDimsLine[7:])

nPntsIdx = [lines.index(l) for l in lines if l.startswith('NPOIN= ')]
nPntsLine = lines[nPntsIdx[0]]
nPnts = int(nPntsLine[7:])

nElemIdx = [lines.index(l) for l in lines if l.startswith('NELEM= ')]
nElemLine = lines[nPntsIdx[0]]
nElem = int(nElemLine[7:])
print("Number of elements: "+ str(nElem) + " \nAs found on line: " + str(nElemIdx[0]))
print(nDims, nPnts)
print(nDimsIdx,nPntsIdx)

arr = np.zeros((nPnts,nDims))

for i in range(nPnts):
    stripLine = lines[i+nPntsIdx[0]+1].strip()
    line = stripLine.split("\t")
    for ii in range(nDims):
        arr[i,ii] = float(line[ii])       
        
x = arr[:,0]
y = arr[:,1]


#x,y = getXY('\\TestMesh.su2')
#xDef,yDef = getXY('/newmesh.su2')


plt.figure()
plt.scatter(x,y)
plt.axis('equal')

#%%

#plt.close('all')
#
#f, (ax1, ax2) = plt.subplots(1, 2)
#ax1.axis('equal')
#
#for i in range(6):    
#    ax1.plot(x[i*6:(i+1)*6],y[i*6:(i+1)*6],color='black')    
#    ax1.plot(x[i:36:6],y[i:36:6], color='black')
#ax1.title.set_text('Initial Mesh')
#
#ax2.axis('equal')
#ax2.xlim([-0.02,1.02])
#for i in range(6):    
#    ax2.plot(xDef[i*6:(i+1)*6],yDef[i*6:(i+1)*6],color='black')    
#    ax2.plot(xDef[i:36:6],yDef[i:36:6], color='black')
#ax2.title.set_text('Deformed Mesh (dx=-0.15, dy=-0.1)')

#plt.figure()
#plt.axis('square')
##plt.scatter(x,y)
#for i in range(6):    
#    plt.plot(x[i*6:(i+1)*6],y[i*6:(i+1)*6],color='black')    
#    plt.plot(x[i:36:6],y[i:36:6], color='black')
#plt.xlim([-0.02,1.02])
#plt.ylim([-0.02,1.02])
#plt.show()


        
