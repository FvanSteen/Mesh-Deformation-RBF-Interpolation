import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import glob

def plot_grid(x,y, ax=None, **kwargs):
    ax = ax or plt.gca()
    segs1 = np.stack((x,y), axis=2)
    segs2 = segs1.transpose(1,0,2)
    ax.add_collection(LineCollection(segs1, **kwargs))
    ax.add_collection(LineCollection(segs2, **kwargs))
    ax.autoscale()

def getXY(fileName):
    os.chdir('../Classical RBF' )
    x = os.getcwd()
    print( "check")
    print(glob.glob(x+ '/*'))
    print(x)
#    fileName = '/newmesh.su2'
    
    fileObj = open(x+fileName, "r") 
    
    
    lines = fileObj.read().splitlines()
    fileObj.close()
    nDimsIdx = [lines.index(l) for l in lines if l.startswith('NDIME= ')]
    nDimsLine = lines[nDimsIdx[0]]
    nDims = int(nDimsLine[7:])
    
    nPntsIdx = [lines.index(l) for l in lines if l.startswith('NPOIN= ')]
    nPntsLine = lines[nPntsIdx[0]]
    nPnts = int(nPntsLine[7:])
    
    
    arr = np.zeros((nPnts,nDims))
    
    for i in range(nPnts):
        line = lines[i+nPntsIdx[0]+1].split("\t")
        for ii in range(nDims):
            arr[i,ii] = float(line[ii])       
            
    x = arr[:,0]
    y = arr[:,1]
    return (x,y)

x,y = getXY('/mesh_NACA0012_inv.su2')
#xDef,yDef = getXY('/newmesh.su2')


plt.figure()
plt.scatter(x,y)


#%%

plt.close('all')

f, (ax1, ax2) = plt.subplots(1, 2)
ax1.axis('equal')

for i in range(6):    
    ax1.plot(x[i*6:(i+1)*6],y[i*6:(i+1)*6],color='black')    
    ax1.plot(x[i:36:6],y[i:36:6], color='black')
ax1.title.set_text('Initial Mesh')

ax2.axis('equal')
ax2.xlim([-0.02,1.02])
for i in range(6):    
    ax2.plot(xDef[i*6:(i+1)*6],yDef[i*6:(i+1)*6],color='black')    
    ax2.plot(xDef[i:36:6],yDef[i:36:6], color='black')
ax2.title.set_text('Deformed Mesh (dx=-0.15, dy=-0.1)')

#plt.figure()
#plt.axis('square')
##plt.scatter(x,y)
#for i in range(6):    
#    plt.plot(x[i*6:(i+1)*6],y[i*6:(i+1)*6],color='black')    
#    plt.plot(x[i:36:6],y[i:36:6], color='black')
#plt.xlim([-0.02,1.02])
#plt.ylim([-0.02,1.02])
#plt.show()


        
