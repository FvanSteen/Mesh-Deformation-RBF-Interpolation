import numpy as np
import matplotlib.pyplot as plt
import os

def polarCartesianTransform(coordsPolar):
    coordsCart = np.empty(np.shape(coordsPolar))
    coordsCart[:,0] = coordsPolar[:,0]*np.cos(coordsPolar[:,1])
    coordsCart[:,1] = coordsPolar[:,0]*np.sin(coordsPolar[:,1])
    return coordsCart    




r_start = 0.2
r_stop = 1.0
theta_start = 0
theta_stop = np.pi/6
z_start = 0
z_stop = 0.2


r = np.linspace(r_start, r_stop, 26)
theta = np.linspace(theta_start, theta_stop, 26)
z = np.linspace(z_start, z_stop, 1)

dims = 2
coords = np.empty((np.size(r)*np.size(theta)*np.size(z),dims))

for k in range(len(z)):
    for j in range(len(theta)):
        for i in range(len(r)):        
            if dims == 2:
                coords[j*(len(r))+i,:] = np.array([r[i],theta[j]])
            elif dims == 3:
                print(k*len(theta)*len(r) + j*len(r) + i)
                coords[k*len(theta)*len(r) + j*len(r) + i,:] = np.array([r[i], theta[j], z[k]])


coordsCart = polarCartesianTransform(coords)

#x_mp = np.array([0.28, 0.44, 0.6, 0.76, 0.92, 0.28, 0.44, 0.6, 0.76, 0.92, 0.199726, 0.199726, 0.199726, 0.199726, 0.199726, 0.99863, 0.99863, 0.99863, 0.99863, 0.99863])
#y_mp = np.array([0, 0, 0, 0, 0, 0.523599, 0.523599, 0.523599, 0.523599, 0.523599, 0.0523599, 0.15708, 0.261799, 0.366519, 0.471239, 0.0523599, 0.15708, 0.261799, 0.366519, 0.471239])
#
#
#idxMoving = np.array([14,15,20,21])
#xdefpol =  np.array([0.0874192 , 0.0874192, 0.0796548,0.0796548 ])
#ydefpol =  np.array([-0.0696985, -0.0696985,-0.0784545, -0.0784545])
#
#
#  
#newPosx= np.array([0.611405,0.770578, 0.604765,0.763695])
#newPosy = np.array([ 0.0951938,  0.118866,   0.184065,   0.211248])
#
#
#
#
#coords[idxMoving,0] = newPosx
#coords[idxMoving,1] = newPosy
#
#coordsCart = polarCartesianTransform(coords)

plt.figure()
plt.scatter(coords[:,0],coords[:,1])
#plt.scatter(coordsCart[:,0],coordsCart[:,1])
plt.scatter(coords[block_idx,0],coords[block_idx,1])
#plt.scatter(coords[17,0],coords[17,1])
#plt.scatter(x_mp,y_mp)
#plt.scatter(coordsCart[idxMoving,0]+xdef, coordsCart[idxMoving,1]+ydef)
#plt.scatter(coords[idxMoving,0], coords[idxMoving,1])

plt.axis('equal')
plt.show()


#%%




#%%



os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\Meshes')
x = os.getcwd()
o = open(x + "\\25x25_per.su2", "w")
#
xElem = 25
yElem = 25
zElem = 0
nBlock = 5;
#
xNodes = xElem +1
yNodes = yElem +1
zNodes = zElem +1
#

block_start_idx = int((yNodes/2-1)*xNodes + xNodes/2-np.floor(nBlock/2)-1)
block_idx = np.arange(block_start_idx,block_start_idx+nBlock+1)
block_idx = np.hstack((block_idx,np.flip(block_idx,axis=0)+xNodes))
block_idx = np.append(block_idx,block_idx[0])

o.write("NDIME= 2 \n")
o.write("NELEM= " + str(xElem*yElem-nBlock) + "\n")

for j in range(yNodes-1):
    for i in range(xNodes-1):
        if(i + j*yNodes in block_idx[0:5]):
            print('skip')
        else:
            o.write("9\t" + str(i + j*yNodes) + "\t" + str(i+1 + j*yNodes) + "\t" +  str(i+1+(1+j)*yNodes) + "\t" +  str(i+(1+j)*yNodes) + "\t" + str(j*xElem+i) +  "\n")
        
o.write("NPOIN= " +  str(xNodes*yNodes*zNodes) + "\n")
for c1,c2 in zip(coordsCart[:,0], coordsCart[:,1]):
    o.write("{0}\t{1}".format(c1,c2) + "\n")
    
o.write("NMARK= 5\n")
o.write("MARKER_TAG= LOWER\n")
o.write("MARKER_ELEMS= "+str(xElem) + "\n")
for i in range(xElem):
    o.write("3\t" + str(i) + "\t" + str(i+1) + "\n")
    
o.write("MARKER_TAG= UPPER\n")
o.write("MARKER_ELEMS= "+str(xElem) + "\n")
for i in range(xElem):
    o.write("3\t" + str(i+yElem*xNodes) + "\t" + str(i+1+yElem*xNodes) + "\n")
    
o.write("MARKER_TAG= LEFT\n")
o.write("MARKER_ELEMS= "+str(yElem) + "\n")
for i in range(yElem):
    o.write("3\t" + str(i*yNodes) + "\t" + str((i+1)*yNodes) + "\n")
    
o.write("MARKER_TAG= RIGHT\n")
o.write("MARKER_ELEMS= "+str(yElem) + "\n")
for i in range(yElem):
    o.write("3\t" + str(i*yNodes+xElem) + "\t" + str((i+1)*yNodes+xElem) + "\n")
    
o.write("MARKER_TAG= BLOCK\n")
o.write("MARKER_ELEMS= "+ str(len(block_idx)-1) + "\n")

for i in range(len( block_idx)-1):
    o.write("3\t" + str(block_idx[i]) + "\t" + str(block_idx[i+1]) + "\n")
#o.write("3\t" + str(block_idx[0]) + "\t" + str(block_idx[-1]) + "\n")
#o.write("3\t" + str(block_idx[-1]) + "\t" + str(block_idx[-1]+xNodes) + "\n")
#o.write("3\t" + str(block_idx[-1]+xNodes) + "\t" + str(block_idx[0]+xNodes) + "\n")
#o.write("3\t" + str(block_idx[0]+xNodes) + "\t" + str(block_idx[0]) + "\n")
#o.write("3\t" + str(14) + "\t" + str(15) + "\n")
#o.write("3\t" + str(15) + "\t" + str(21) + "\n")
#o.write("3\t" + str(21) + "\t" + str(20) + "\n")
#o.write("3\t" + str(20) + "\t" + str(14) + "\n")
o.close()

#%%
from meshpy.tet import MeshInfo, build

mesh_info = MeshInfo()
mesh_info.set_points([
    (0,0,0), (2,0,0), (2,2,0), (0,2,0),
    (0,0,12), (2,0,12), (2,2,12), (0,2,12),
    ])
mesh_info.set_facets([
    [0,1,2,3],
    [4,5,6,7],
    [0,4,5,1],
    [1,5,6,2],
    [2,6,7,3],
    [3,7,4,0],
    ])
mesh = build(mesh_info)
print("Mesh Points:")
for i, p in enumerate(mesh.points):
    print( i, p)
print("Point numbers in tetrahedra:")
for i, t in enumerate(mesh.elements):
    print(i, t)