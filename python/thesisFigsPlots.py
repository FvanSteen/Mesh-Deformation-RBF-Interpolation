# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 11:16:22 2023

@author: floyd
"""


plt.close('all')
#x = [0,1,0.5,0]
#y = [0,0,0.5,0]
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
a = plt.plot(x[0:-1],y[0:-1], color = 'black', label = 'Boundary line segments')
b = plt.scatter(x[:-1],y[:-1],color = 'black', label = 'Sliding edge nodes')
c = plt.scatter(midpntx[:-1],midpnty[:-1],color = 'red', zorder=3, label = 'Line midpoints')
for i in range(len(n_x)-1):
#i = 
    if i == 0:
        arrow1 = plt.arrow(midpntx[i], midpnty[i],n_x[i],n_y[i], head_width = .02, color='red', zorder = 3, label = 'Line normal vector')
    else:
        plt.arrow(midpntx[i], midpnty[i],n_x[i],n_y[i], head_width = .02, color='red', zorder = 3, label = '_nolegend_')
#    plt.arrow([midpntx[i], midpntx[i]+n_x[i]],[midpnty[i], midpnty[i]+n_y[i]], color='blue')
#    plt.plot([x[i], x[i]+t_x[i]],[y[i], y[i]+t_y[i]], color='blue')
d = plt.scatter(dx[2],dy[2], color = 'blue', zorder=3, label = 'Displaced sliding node')    

xn = dx+5*n_x[:-1]*(deltax*(n_x[:-1]/np.sqrt(n_x[:-1]**2+n_y[:-1]**2)) + deltay*(n_y[:-1]/np.sqrt(n_x[:-1]**2+n_y[:-1]**2)))
yn = dy+5*n_y[:-1]*(deltax*(n_x[:-1]/np.sqrt(n_x[:-1]**2+n_y[:-1]**2)) + deltay*(n_y[:-1]/np.sqrt(n_x[:-1]**2+n_y[:-1]**2)))
#for i in range(len(deltax)):
i =2
e = plt.arrow(dx[i],dy[i], 3*n_x[i]*(deltax[i]*(n_x[i]/np.sqrt(n_x[i]**2+n_y[i]**2)) + deltay[i]*(n_y[i]/np.sqrt(n_x[i]**2+n_y[i]**2))), 3*n_y[i]*(deltax[i]*(n_x[i]/np.sqrt(n_x[i]**2+n_y[i]**2)) + deltay[i]*(n_y[i]/np.sqrt(n_x[i]**2+n_y[i]**2))),  head_width = .02, color = 'blue', label = 'Projection')
#plt.scatter(xn,yn, color = 'green', zorder = 3)
f = plt.scatter(xn[2],yn[2],color = 'blue', facecolor = 'none')
plt.axis('equal')
plt.ylim((-0.1,1.25))
plt.axis('off')

plt.legend([a,b,c,arrow1,d,e, f], ['Boundary line segments','Sliding edge nodes','Line midpoints','Line normal vector','Displaced sliding node', 'Projection', 'Projected node'], handler_map={mpatches.FancyArrow : HandlerPatch(patch_func=make_legend_arrow),}, loc = 'lower right')
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
a = plt.plot(x,y, color = 'black')
b = plt.scatter(x,y, color = 'black')
c = plt.scatter(mpx,mpy, color = 'red',zorder = 3)

for i in range(len(n_x)):
#i = 
    if i == 0:
        arrow1 = plt.arrow(mpx[i], mpy[i],-n_x[i],-n_y[i], head_width = .07, color='red', zorder = 3, label = 'Line normal vector')
    else:
        plt.arrow(mpx[i], mpy[i],-n_x[i],-n_y[i], head_width = .07, color='red', zorder = 3, label = '_nolegend_')
        
#d = plt.scatter(-1.25,1.3,color = 'blue')
#e = plt.arrow(-1.25,1.3,1.1,0, head_width = .07, color='blue', zorder = 3)
#f = plt.scatter(0,1.3,color = 'blue', marker = 'o',facecolor='none')

d = plt.scatter(0,1.3,color = 'blue')
e = plt.arrow(0,1.3,0.05,-0.05, head_width = .07, color='blue', zorder = 3)
f = plt.scatter(0.15,1.3-0.15,color = 'blue', marker = 'o',facecolor='none')
plt.legend([a,b,c,arrow1,d,e,f], ['Boundary line segments','Sliding edge nodes','Line midpoints','Line normal vector','Displaced sliding node', 'Projection', 'Projected node'], handler_map={mpatches.FancyArrow : HandlerPatch(patch_func=make_legend_arrow),}, loc = 'lower right')
plt.axis('equal')
plt.ylim((-0.5,3.5))
plt.xlim((-1.5,3))
plt.axis('off')
plt.show()
plt.savefig(figPath + "projection3.png", dpi=800,bbox_inches='tight')