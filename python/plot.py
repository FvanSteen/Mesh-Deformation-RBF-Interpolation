import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

def plot_grid(x,y, ax=None, **kwargs):
    ax = ax or plt.gca()
    segs1 = np.stack((x,y), axis=2)
    segs2 = segs1.transpose(1,0,2)
    ax.add_collection(LineCollection(segs1, **kwargs))
    ax.add_collection(LineCollection(segs2, **kwargs))
    ax.autoscale()

os.chdir('../Classical RBF' )
x = os.getcwd()

fileName = '/newmesh.su2'

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


plt.close('all')

plt.scatter(x,y)


#%%
fig = plt.figure()
ax = fig.add_subplot(111)
if nDims ==2:
    ax.scatter(arr[:,0],arr[:,1])
elif nDims ==3:
    ax.scatter(arr[:,0],arr[:,1],arr[:,2])    
plt.show()
        
