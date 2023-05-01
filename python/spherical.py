import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import os

def cylindricalCartesianTransform(coordsCylindrical):
    coordsCart = np.empty(np.shape(coordsCylindrical))
    coordsCart[:,0] = coordsCylindrical[:,0]*np.cos(coordsCylindrical[:,1])
    coordsCart[:,1] = coordsCylindrical[:,0]*np.sin(coordsCylindrical[:,1])
    coordsCart[:,2] = coordsCylindrical[:,2]
    return coordsCart    




r_start = 0.2
r_stop = 1.0
theta_start = 0
theta_stop = np.pi/6
z_start = 0
z_stop = 0.5


r = np.linspace(r_start, r_stop, 6)
theta = np.linspace(theta_start, theta_stop, 6)
z = np.linspace(z_start, z_stop, 6)

dims = 3
coords = np.empty((np.size(r)*np.size(theta)*np.size(z),dims))

for k in range(len(z)):
    for j in range(len(theta)):
        for i in range(len(r)):        
            if dims == 2:
                coords[j*(len(r))+i,:] = np.array([r[i],theta[j]])
            elif dims == 3:
                coords[k*len(theta)*len(r) + j*len(r) + i,:] = np.array([r[i], theta[j], z[k]])

coordsCart = cylindricalCartesianTransform(coords)

#%% CARTESIAN COORDINATES
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(coordsCart[:,0],coordsCart[:,1],coordsCart[:,2],marker='1')
scaling = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
ax.auto_scale_xyz(*[[np.min(scaling), np.max(scaling)]]*3)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()
#%% CYLINDRICAL COORDINATES

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(coords[:,0],coords[:,1],coords[:,2],marker='1')
scaling = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
ax.auto_scale_xyz(*[[np.min(scaling), np.max(scaling)]]*3)
ax.set_xlabel('r')
ax.set_ylabel('theta')
ax.set_zlabel('z')
plt.show()



#%%


os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\Meshes')
x = os.getcwd()
o = open(x + "\\5x5x5_per.su2", "w")
#
xElem = np.size(r)-1
yElem = np.size(theta)-1
zElem = np.size(z)-1
#
xNodes = xElem +1
yNodes = yElem +1
zNodes = zElem +1

middleElem = int((xElem*yElem*zElem-1)/2)
startingNode = int((zNodes/2-1)*xNodes*yNodes + (yNodes/2-1)*xNodes + xNodes/2-1)
#
o.write("NDIME= 3 \n")
o.write("NELEM= " + str(xElem*yElem*zElem-1) + "\n")
for k in range(zNodes-1):
    for j in range(yNodes-1):
        for i in range(xNodes-1):
            if(k*xElem*yElem+ j*xElem+i != middleElem): 
#                o.write("12\t" + str(i + j*yNodes + k*yNodes*zNodes) + "\t" + str(i+1 + j*yNodes + k*yNodes*zNodes) + "\t" +  str(i+1+(1+j)*yNodes) + "\t" +  str(i+(1+j)*yNodes + k*yNodes*zNodes) + "\t" + str(j*xElem+i + k*yNodes*zNodes) +  "\n")
#                print(i + j*yNodes + k*yNodes*xNodes, i+1 + j*yNodes + k*yNodes*xNodes, i+1 + (j+1)*yNodes + k*yNodes*xNodes, i + (j+1)*yNodes + k*yNodes*xNodes, i + j*yNodes + (k+1)*yNodes*xNodes, i+1 + j*yNodes + (k+1)*yNodes*xNodes,  i+1 + (j+1)*yNodes + (k+1)*yNodes*xNodes, i + (j+1)*yNodes + (k+1)*yNodes*xNodes)
                print("12\t" + str(i + j*yNodes + k*yNodes*xNodes) + "\t" +  str(i+1 + j*yNodes + k*yNodes*xNodes)  + "\t" + str(i+1 + (j+1)*yNodes + k*yNodes*xNodes)  + "\t" + str(i + (j+1)*yNodes + k*yNodes*xNodes)  + "\t" +  str(i + j*yNodes + (k+1)*yNodes*xNodes)  + "\t" +  str(i+1 + j*yNodes + (k+1)*yNodes*xNodes)  + "\t" +  str(i+1 + (j+1)*yNodes + (k+1)*yNodes*xNodes)  + "\t" +  str(i + (j+1)*yNodes + (k+1)*yNodes*xNodes) + "\n")
                o.write("12\t" + str(i + j*yNodes + k*yNodes*xNodes) + "\t" +  str(i+1 + j*yNodes + k*yNodes*xNodes)  + "\t" + str(i+1 + (j+1)*yNodes + k*yNodes*xNodes)  + "\t" + str(i + (j+1)*yNodes + k*yNodes*xNodes)  + "\t" +  str(i + j*yNodes + (k+1)*yNodes*xNodes)  + "\t" +  str(i+1 + j*yNodes + (k+1)*yNodes*xNodes)  + "\t" +  str(i+1 + (j+1)*yNodes + (k+1)*yNodes*xNodes)  + "\t" +  str(i + (j+1)*yNodes + (k+1)*yNodes*xNodes) + "\n")

            
o.write("NPOIN= " +  str(xNodes*yNodes*zNodes) + "\n")
for c1,c2,c3 in zip(coordsCart[:,0], coordsCart[:,1], coordsCart[:,2]):
    o.write("{0}\t{1}\t{2}".format(c1,c2,c3) + "\n")
    
o.write("NMARK= 7\n")
o.write("MARKER_TAG= LOWER\n")
o.write("MARKER_ELEMS= "+str(xElem*yElem) + "\n")
for j in range(yElem):
    for i in range(xElem):
        o.write("9\t" + str(i+j*xNodes) + "\t" + str(i+j*xNodes+1) + "\t" + str(i+(j+1)*xNodes+1) + "\t" + str(i+(j+1)*xNodes) + "\n")
    
o.write("MARKER_TAG= UPPER\n")
o.write("MARKER_ELEMS= "+str(xElem*yElem) + "\n")
for j in range(yElem):
    for i in range(xElem):    
        o.write("9\t" + str(i+j*xNodes + (zElem)*xNodes*yNodes) + "\t" + str(i+j*xNodes+1+ (zElem)*xNodes*yNodes) + "\t" + str(i+(j+1)*xNodes+1+ (zElem)*xNodes*yNodes) + "\t" + str(i+(j+1)*xNodes+ (zElem)*xNodes*yNodes) + "\n")
    
o.write("MARKER_TAG= LEFT\n")
o.write("MARKER_ELEMS= "+str(yElem*zElem) + "\n")
for j in range(zElem):
    for i in range(yElem):
        o.write("9\t" + str(i*xNodes + j*xNodes*yNodes) + "\t" + str((i+1)*xNodes+ j*xNodes*yNodes) + "\t" + str((i+1)*xNodes+ (j+1)*xNodes*yNodes) + "\t" + str((i)*xNodes+ (j+1)*xNodes*yNodes)+ "\n")

o.write("MARKER_TAG= RIGHT\n") 
o.write("MARKER_ELEMS= "+str(yElem*zElem) + "\n")
for j in range(zElem):
    for i in range(yElem):
        o.write("9\t" + str(i*xNodes + j*xNodes*yNodes + xElem) + "\t" + str((i+1)*xNodes+ j*xNodes*yNodes+ xElem) + "\t" + str((i+1)*xNodes+ (j+1)*xNodes*yNodes+ xElem) + "\t" + str((i)*xNodes+ (j+1)*xNodes*yNodes+ xElem)+ "\n")
    
o.write("MARKER_TAG= FRONT\n") 
o.write("MARKER_ELEMS= "+str(xElem*zElem) + "\n")
for j in range(zElem):
    for i in range(xElem):
        o.write("9\t" + str(i + j*xNodes*yNodes) + "\t" + str(i+1 + j*xNodes*yNodes) + "\t" + str(i+1 + (j+1)*xNodes*yNodes ) + "\t" + str(i + (j+1)*xNodes*yNodes)+ "\n")
        
o.write("MARKER_TAG= BACK\n") 
o.write("MARKER_ELEMS= "+str(xElem*zElem) + "\n")
for j in range(zElem):
    for i in range(xElem):
        o.write("9\t" + str(i + j*xNodes*yNodes + xNodes*yElem) + "\t" + str(i+1 + j*xNodes*yNodes+ xNodes*yElem) + "\t" + str(i+1 + (j+1)*xNodes*yNodes + xNodes*yElem) + "\t" + str(i + (j+1)*xNodes*yNodes+ xNodes*yElem)+ "\n")
    

    
o.write("MARKER_TAG= BLOCK\n")
o.write("MARKER_ELEMS= "+ str(6) + "\n")
o.write("9\t" + str(startingNode) + "\t" + str(startingNode+1) + "\t" + str(startingNode+1+xNodes) + "\t" + str(startingNode+xNodes)  + "\n")
o.write("9\t" + str(startingNode+xNodes*yNodes) + "\t" + str(startingNode+1+xNodes*yNodes) + "\t" + str(startingNode+1+xNodes+xNodes*yNodes) + "\t" + str(startingNode+xNodes+xNodes*yNodes)  + "\n")
o.write("9\t" + str(startingNode) + "\t" + str(startingNode+xNodes*yNodes) + "\t" + str(startingNode+xNodes*yNodes+xNodes) + "\t" + str(startingNode+xNodes)  + "\n")
o.write("9\t" + str(startingNode+1) + "\t" + str(startingNode+xNodes*yNodes+1) + "\t" + str(startingNode+xNodes*yNodes+xNodes+1) + "\t" + str(startingNode+xNodes+1)  + "\n")
o.write("9\t" + str(startingNode) + "\t" + str(startingNode+1) + "\t" + str(startingNode+1+xNodes*yNodes) + "\t" + str(startingNode+xNodes*yNodes)  + "\n")
o.write("9\t" + str(startingNode+xNodes) + "\t" + str(startingNode+1+xNodes) + "\t" + str(startingNode+1+xNodes*yNodes+xNodes) + "\t" + str(startingNode+xNodes*yNodes+xNodes)  + "\n")
o.close()