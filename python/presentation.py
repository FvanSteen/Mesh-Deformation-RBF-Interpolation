# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 17:31:29 2023

@author: floyd
"""

#import numpy as np
#import matplotlib.pyplot as plt
#
#nodes = np.linspace(-2,2,9)
#zero = np.zeros((9,1))
#
#alpha = 0.25
#r = 2
#xi = abs(nodes/r)
#f = (1-xi)**4*(4*xi+1)
#f[xi>1] = 0
#
#delta = alpha*f
#
#nodes2 = np.linspace(-2,2,101)
#xi2 = abs(nodes2/r)
#f = (1-xi2)**4*(4*xi2+1)
#
#plt.close("all")
#
#
#ax1 = plt.subplot(212)
#plt.scatter(nodes,zero,label='Initial')
##plt.ylim([-0.005,0.005])
#
#plt.scatter(nodes+delta,zero+0.01, label='deformed')
#plt.scatter([0,0.25],[0,0.01],label = 'Control node',color = 'red')
#ax1.axes.get_yaxis().set_visible(False)
#plt.legend()
#
##plt.ylim([-0.005,0.005])
#ax2 = plt.subplot(211, sharex=ax1)
#plt.plot(nodes2,f)
#plt.title("Wendland $C^2$, $r=2$")
#plt.show()


#import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib.animation as animation
#from colMap import colMap
#from mpl_toolkits.mplot3d.art3d import Poly3DCollection
#from funs import getMeshQual, getPlotData
#import matplotlib.collections
#from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
#import os
#os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\Meshes')
#cmapMatlab = colMap()
#fig = plt.figure()
#
#
#
#ims = []
#fileName = '/25x25.su2' 
#[f,v,elemType,bdryPnts,markerTags, nElemsMarks, FFD_pnts] = getPlotData(fileName) 
#for i in range(20):
#    meshQual = getMeshQual(fileName)       
#    quad_idx = np.where(elemType == 9)
#    colors = cmapMatlab(plt.Normalize(0,1)(meshQual))
#    pc_quad = matplotlib.collections.PolyCollection(v[f[quad_idx][:,0:4]],cmap=cmapMatlab,  facecolors=colors, edgecolor="black",linewidth=0.1)
#    
#    ax = fig.add_subplot(1,1,1)
#    
#    
#    polys = ax.add_collection(pc_quad)
#    
#    ax.autoscale()
#    ax.set_aspect('equal')
# 
#    
#    polys.set_clim(0,1)
#    
#    aspect = 20
#    pad_fraction = 0.5
#    divider = make_axes_locatable(ax)
#    width = axes_size.AxesY(ax, aspect=1./aspect)
#    pad = axes_size.Fraction(pad_fraction, width)
#    cax = divider.append_axes("right", size=width, pad=pad)
#    plt.colorbar(polys, cax=cax)
#    
#    
#    im = plt.show()#plt.imshow(f(x, y), animated=True)
#    ims.append([im])
#
#ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True,
#                                repeat_delay=1000)
#
#ani.save('dynamic_images.gif')
#
#plt.show()     

#%%


import imageio
pdir = "C:/Users/floyd/git/Mesh-Deformation-RBF-Interpolation/python/figs"
filenames = ["/sq_none_0","/sq_none_1","/sq_none_2","/sq_none_3","/sq_none_4","/sq_none_5","/sq_none_6","/sq_none_7","/sq_none_8","/sq_none_9"]
images = []
#for filename in filenames:
for i in range(20):
    filename = "/sq_ds_"+str(i) 
    images.append(imageio.imread(pdir+filename+".png"))
    
for _ in range(10):
    images.append(imageio.imread(pdir+filename+".png"))
imageio.mimsave('C:/Users/floyd/git/Mesh-Deformation-RBF-Interpolation/python/ds.gif', images)

