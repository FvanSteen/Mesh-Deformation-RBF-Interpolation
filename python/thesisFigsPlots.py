# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 11:16:22 2023

@author: floyd
"""
import os
import numpy as np
from matplotlib import pyplot as plt
figPath = os.path.dirname(os.path.abspath(__file__)) + "/figs/"
plt.close('all')
#x = [0,1,0.5,0]
#y = [0,0,0.5,0]

#%%
plt.close('all')

x = np.arange(0,5,0.01)

phi = (1-abs(x))**4 * (4*abs(x)+1)
phi_r2 = (1-abs(x/2))**4 * (4*abs(x/2)+1)
phi_r5 = (1-abs(x/5))**4 * (4*abs(x/5)+1)
phi[abs(x)>1] = 0
phi_r2[abs(x)>2] = 0
phi_r5[abs(x)>5] = 0

plt.figure()
plt.plot(x,phi,linewidth=1.5)
plt.plot(x,phi_r2,linewidth=1.5)
plt.plot(x,phi_r5,linewidth=1.5)
plt.xlabel('$\\xi$ $[-]$',size=12)
plt.ylabel('$\phi(\\xi)$ $[-]$',size=12)
plt.legend(['$r=1$','$r=2$','$r=5$'], fontsize =12)
plt.savefig(figPath + "rbf_xi.png", dpi=800,bbox_inches = 'tight')
#%%
from matplotlib import cm
from colMap import colMap
X = np.arange(-1.5, 1.5, 0.01)
Y = np.arange(-1.5, 1.5, 0.01)
X, Y = np.meshgrid(X, Y)
r = np.sqrt(X**2 + Y**2)
support = 1


Z = (1-r)**4 * (4*r+1)
Z[r/support>1] = 0
cmapMatlab = colMap();
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(X, Y, Z, cmap=cm.coolwarm, vmin=0,vmax=1, linewidth=0, antialiased=False)
plt.xlabel('$x$',size=12)
plt.ylabel('$y$', size=12)
ax.zaxis.labelpad=10
ax.zaxis.set_rotate_label(False)
ax.set_zlabel('$\phi(\\xi)$', size =12, rotation = 0)
v = np.linspace(0, 1.0, 6, endpoint=True)
fig.colorbar(surf, shrink=0.5, ticks = v)
ax.view_init(elev=30., azim=65)
plt.savefig(figPath + "rbf_xy-plane.png", dpi=800)


#%%
x = np.linspace(0,1,100)
l = 1
y = l/np.pi*np.sin(x*np.pi/l)
plt.figure()
plt.plot(x,y, linewidth=1.5)


plt.axis([0,1,0,1])
plt.gca().set_aspect('equal', adjustable='box')
plt.xlabel('$y-y_c$ $[m]$ ', size = 12)
plt.ylabel('$\delta_y$ $[m]$', size = 12, )
plt.show()

plt.savefig(figPath + "periodic_sine.png", dpi=800,bbox_inches='tight')


#%%
x = np.array([[0,0],[1,0],[0.5,.5*np.sqrt(2)],[0,0]])

plt.close('all')
plt.figure()

for i in range(np.shape(x)[0]-1):
    plt.plot( x[i:i+2,0], x[i:i+2,1], color= 'black')

delta = 0.03
plt.text(-1.5*delta,-delta,'0',fontsize = 14)
plt.text(x[1,0]+delta,x[1,1]-delta,'1', fontsize = 14)
plt.text(x[2,0]+delta/2,x[2,1]+delta/2,'2', fontsize = 14)
plt.axis('equal')
plt.axis('off')
plt.show()
plt.savefig(figPath + "triangle.png", dpi=800,bbox_inches='tight')



#%%
x = np.array([[0,0],[1,0],[1,1],[0,1],[0,0]])
plt.figure()
for i in range(np.shape(x)[0]-1):
    plt.plot( x[i:i+2,0], x[i:i+2,1], color= 'black')

delta = 0.03
plt.text(-2*delta,-delta,'0',fontsize = 14)
plt.text(x[1,0]+delta,x[1,1]-delta,'1', fontsize = 14)
plt.text(x[2,0]+delta/2,x[2,1]+delta/2,'2', fontsize = 14)
plt.text(x[3,0]-delta,x[3,1]+delta/2,'3', fontsize = 14)
plt.axis('equal')
plt.axis('off')
plt.show()
plt.savefig(figPath + "Quad.png", dpi=800,bbox_inches='tight')


#%%

x = np.array([[0,0],[0.9,-0.15],[0.5,0.8*.5*np.sqrt(2)],[0,0]])
x2 = np.array([[0.9,-0.15],[0.95,0.17],[0.5,0.8*.5*np.sqrt(2)],[0,0]])
x3 = np.array([[0,0],[0.95,0.17]])

plt.close('all')
plt.figure()

for i in range(np.shape(x)[0]-1):
    plt.plot( x[i:i+2,0], x[i:i+2,1], color= 'black')
    plt.plot( x2[i:i+2,0], x2[i:i+2,1], color= 'black')
    plt.plot( x3[i:i+2,0], x3[i:i+2,1], color= 'black', linestyle = '--')

delta = 0.03
plt.text(-1.5*delta,-delta,'0',fontsize = 14)
plt.text(x[1,0]+delta,x[1,1]-delta,'1', fontsize = 14)
plt.text(x[2,0]+delta/2,x[2,1]+delta/2,'3', fontsize = 14)
plt.text(x2[1,0]+delta/2,x2[1,1]+delta/2,'2', fontsize = 14)
plt.axis('equal')
plt.axis('off')
plt.show()
plt.savefig(figPath + "Tetrahedral.png", dpi=800,bbox_inches='tight')
#%%
plt.close('all')
x = np.array([[0,0],[1,0],[1,1],[0,1],[0,0]])

x2 = x+0.15
x2[:,0]  = x2[:,0] + 0.05
plt.figure()
for i in range(np.shape(x)[0]-1):
    plt.plot( x[i:i+2,0], x[i:i+2,1], color= 'black')
#    plt.plot( x2[i:i+2,0], x2[i:i+2,1], color= 'black')
plt.plot(x2[1:4,0], x2[1:4,1], color= 'black')
plt.plot( [x[0,0], x2[0,0]], [x[0,1], x2[0,1]] , color= 'black', linestyle = '--')
plt.plot( [x2[3,0], x2[0,0]], [x2[3,1], x2[0,1]] , color= 'black', linestyle = '--')
plt.plot( [x2[1,0], x2[0,0]], [x2[1,1], x2[0,1]] , color= 'black', linestyle = '--')
plt.plot( [x[1,0], x2[1,0]], [x[1,1], x2[1,1]] , color= 'black')
plt.plot( [x[2,0], x2[2,0]], [x[2,1], x2[2,1]] , color= 'black')
plt.plot( [x[3,0], x2[3,0]], [x[3,1], x2[3,1]] , color= 'black')
delta = 0.03
plt.text(-2*delta,-delta,'0',fontsize = 14)
plt.text(x[1,0]+delta,x[1,1]-1.5*delta,'1', fontsize = 14)
plt.text(x[2,0]-delta,x[2,1]+delta,'5', fontsize = 14)
plt.text(x[3,0]-2*delta,x[3,1]+delta/2,'4', fontsize = 14)

plt.text(x2[0,0]+0.75*delta,x2[0,1]+0.75*delta,'3', fontsize = 14)
plt.text(x2[1,0]+delta,x2[1,1]+0.5*delta,'2', fontsize = 14)
plt.text(x2[2,0]+0.5*delta,x2[2,1]+0.5*delta,'6', fontsize = 14)
plt.text(x2[3,0]-1.5*delta,x2[3,1]+delta/2,'7', fontsize = 14)
plt.axis('equal')
plt.axis('off')
plt.show()
plt.savefig(figPath + "hexahedral.png", dpi=800,bbox_inches='tight')

#%%
from matplotlib.legend_handler import HandlerPatch
import matplotlib.patches as mpatches

def make_legend_arrow(legend, orig_handle,
                      xdescent, ydescent,
                      width, height, fontsize):
    p = mpatches.FancyArrow(0, 0.5*height, width, 0, length_includes_head=True, head_width=0.75*height )
    return p

    
r = 1;
theta = np.linspace(90,202.5,6)

x = r*np.cos(theta/180*np.pi)
y = r*np.sin(theta/180*np.pi)

midpntx = (x[1:]+x[0:-1])/2
midpnty = (y[1:]+y[0:-1])/2
    
t_x = x[1:]-x[0:-1]
t_y = y[1:]-y[0:-1]
n_x = 0.5*t_y
n_y = 0.5*-t_x

dx = x[1:-1] - .75*t_x[1:]
dy = y[1:-1] - .75*t_y[1:]

deltax = (midpntx[:-1]-dx)
deltay = (midpnty[:-1]-dy)

plt.figure()
a = plt.plot(x[0:-1],y[0:-1], color = 'black', label = 'Boundary line segments',linewidth = .8)
b = plt.scatter(x[:-1],y[:-1],color = 'black', label = 'Sliding edge nodes', s=20)
c = plt.scatter(midpntx[2],midpnty[2],color = 'C1', zorder=3, label = 'Edge midpoints')
for i in range(len(n_x)-1):
#i = 
    if i == 2:
        arrow1 = plt.arrow(midpntx[i], midpnty[i],n_x[i],n_y[i], head_width = .02, color='C1', zorder = 3, label = 'Edge normal vector')
#    else:
#        plt.arrow(midpntx[i], midpnty[i],n_x[i],n_y[i], head_width = .02, color='C1', zorder = 3, label = '_nolegend_')
#    plt.arrow([midpntx[i], midpntx[i]+n_x[i]],[midpnty[i], midpnty[i]+n_y[i]], color='blue')
#    plt.plot([x[i], x[i]+t_x[i]],[y[i], y[i]+t_y[i]], color='blue')
d = plt.scatter(dx[2],dy[2], color = 'C0', zorder=3, label = 'Displaced sliding node')    

xn = dx+5*n_x[:-1]*(deltax*(n_x[:-1]/np.sqrt(n_x[:-1]**2+n_y[:-1]**2)) + deltay*(n_y[:-1]/np.sqrt(n_x[:-1]**2+n_y[:-1]**2)))
yn = dy+5*n_y[:-1]*(deltax*(n_x[:-1]/np.sqrt(n_x[:-1]**2+n_y[:-1]**2)) + deltay*(n_y[:-1]/np.sqrt(n_x[:-1]**2+n_y[:-1]**2)))
#for i in range(len(deltax)):
i =2
e = plt.arrow(dx[i],dy[i], 3*n_x[i]*(deltax[i]*(n_x[i]/np.sqrt(n_x[i]**2+n_y[i]**2)) + deltay[i]*(n_y[i]/np.sqrt(n_x[i]**2+n_y[i]**2))), 3*n_y[i]*(deltax[i]*(n_x[i]/np.sqrt(n_x[i]**2+n_y[i]**2)) + deltay[i]*(n_y[i]/np.sqrt(n_x[i]**2+n_y[i]**2))),  head_width = .02, color = 'C0', label = 'Projection')
#plt.scatter(xn,yn, color = 'green', zorder = 3)
f = plt.scatter(xn[2],yn[2],color = 'C0', facecolor = 'none', zorder = 4)
plt.axis('equal')
plt.ylim((-0.1,1.25))
plt.axis('off')

plt.legend([a,b,c,arrow1,d,e, f], ['Boundary edge segments','Sliding edge nodes','Nearest edge midpoint','Edge normal vector','Displaced sliding node', 'Projection', 'Projected node'], handler_map={mpatches.FancyArrow : HandlerPatch(patch_func=make_legend_arrow),}, loc = 'lower right')
plt.show()
plt.savefig(figPath + "projection.png", dpi=800,bbox_inches='tight')


#%%
plt.close('all')
x = np.array([0,0,0.8*np.sqrt(2),1+0.8*np.sqrt(2)])
y = np.array([0,1,1+0.8*np.sqrt(2), 1+0.8*np.sqrt(2)])

mpx = (x[1:]+x[:-1])/2
mpy = (y[1:]+y[:-1])/2

t_x = x[1:]-x[0:-1]
t_y = y[1:]-y[0:-1]
n_x = 0.5*t_y
n_y = 0.5*-t_x

plt.figure()
a = plt.plot(x,y, color = 'black', linewidth=.8)
b = plt.scatter(x,y, color = 'black', s = 20)
c = plt.scatter(mpx[0],mpy[0], color = 'C1',zorder = 3)

for i in range(len(n_x)):
#i = 
    if i == 0:
        arrow1 = plt.arrow(mpx[i], mpy[i],-n_x[i],-n_y[i], head_width = .07, color='C1', zorder = 3, label = 'Line normal vector')
#    else:
#        plt.arrow(mpx[i], mpy[i],-n_x[i],-n_y[i], head_width = .07, color='C1', zorder = 3, label = '_nolegend_')
        
d = plt.scatter(-1.25,1.3,color = 'C0')
e = plt.arrow(-1.25,1.3,1.1,0, head_width = .07, color='C0', zorder = 3)
f = plt.scatter(0,1.3,color = 'C0', marker = 'o',facecolor='none')

#d = plt.scatter(0,1.3,color = 'C0')
#e = plt.arrow(0,1.3,0.05,-0.05, head_width = .07, color='C0', zorder = 3)
#f = plt.scatter(0.15,1.3-0.15,color = 'C0', marker = 'o',facecolor='none')
plt.legend([a,b,c,arrow1,d,e,f], ['Boundary edge segments','Sliding edge nodes','Nearest edge midpoint','Edge normal vector','Displaced sliding node', 'Projection', 'Projected node'], handler_map={mpatches.FancyArrow : HandlerPatch(patch_func=make_legend_arrow),}, loc = 'lower right')
plt.axis('equal')
plt.ylim((-0.5,3.5))
plt.xlim((-1.5,3))
plt.axis('off')
plt.show()
plt.savefig(figPath + "projection2.png", dpi=800,bbox_inches='tight')