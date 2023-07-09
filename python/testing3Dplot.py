# -*- coding: utf-8 -*-
"""
Created on Sun Oct 16 14:28:47 2022

@author: floyd
"""

# Import libraries
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
 
 
# Create axis
axes = [5, 5, 5]
 
# Create Data
data = np.ones(axes, dtype=np.bool)
 
# Control Transparency
alpha = 0.9
 
# Control colour
colors = np.empty(axes + [4], dtype=np.float32)
 

colors[0] = [1, 0, 0, alpha]  # red
colors[1] = [0, 1, 0, alpha]  # green
colors[2] = [0, 0, 1, alpha]  # blue
colors[3] = [1, 1, 0, alpha]  # yellow
colors[4] = [1, 1, 1, alpha]  # grey
# Plot figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
 
# Voxels is used to customizations of the
# sizes, positions and colors.
ax.voxels(data, facecolors=colors)

#%%

from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
fig = plt.figure()
ax = Axes3D(fig)
x = [0,1,1,0]
y = [0,0,1,1]
z = [0,1,0,1]
verts = [list(zip(x,y,z))]
ax.add_collection3d(Poly3DCollection(verts, edgecolor='black'))
plt.show()
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

#%%
#%%
#import numpy as np
#tarr = []
##tarr.reshape((1,3))
#
#tarr.append(np.transpose(np.array([[0,0,0], [1,1,1]])))
##tarr.append(np.array([[0,0,0], [1,1,1]]))
#print(np.array(tarr))


#%%

# Coordinates of all points according to the order given in the mesh file
#v = np.array([[0,0,0],
#              [1,0,0],
#              [2,0,0],
#              [0,1,0],
#              [1,1,0],
#              [2,1,0],
#              [0,0,1],
#              [1,0,1],
#              [2,0,1],
#              [0,1,1],
#              [1,1,1],
#              [2,1,1]])
#
## f contains the indices that correspond to an element
#f = np.array([[0,1,7,6,3,4,10,9],
#              [1,2,8,7,4,5,11,10]])
#        
#i = 0
#verts = [np.array([v[f][i][0],v[f][i][1],v[f][i][5],v[f][i][4]]),
#          np.array([v[f][i][3],v[f][i][2],v[f][i][6],v[f][i][7]]),
#          np.array([v[f][i][3],v[f][i][0],v[f][i][4],v[f][i][7]]),
#          np.array([v[f][i][2],v[f][i][1],v[f][i][5],v[f][i][6]]),
#          np.array([v[f][i][0],v[f][i][1],v[f][i][2],v[f][i][3]]),
#          np.array([v[f][i][4],v[f][i][5],v[f][i][6],v[f][i][7]])]
#
#
#plt.close('all')
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1,projection='3d')
##for i in range(len(f)):
#for i in range():
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

