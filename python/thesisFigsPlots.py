import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.patches import Wedge
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import matplotlib.lines as mlines

figPath = os.path.dirname(os.path.abspath(__file__)) + "/figs/"
#%%
plt.figure()
plt.rcParams["figure.figsize"] = [6.4, 4.8]
plt.plot([0,1,1,0,0],[0,0,1,1,0], linewidth = 2, color = 'black',label = '_nolegend_')
plt.plot([1,2,2,1,1],[0,0,1,1,0], linewidth = 2, color = 'black',label = '_nolegend_')
plt.plot([0,1],[0.5,0.5],color = 'black',linewidth=2,label = '_nolegend_')
plt.scatter([0,1,1,0,0,2,2],[0,0,1,1,0.5,0,1], color = 'black', s = 75,label = '_nolegend_')
plt.scatter(1,0.5, color = 'red',s =75, zorder = 4)
plt.legend(['Hanging node'], bbox_to_anchor=(0.975,0.94), loc = 'upper right')
plt.axis('equal')
plt.axis('off')
plt.savefig(figPath + "hanging_node.png", dpi=150 ,bbox_inches='tight')

#%%
plt.close("all")
from matplotlib.legend_handler import HandlerPatch

def make_legend_arrow(legend, orig_handle,
                      xdescent, ydescent,
                      width, height, fontsize):
    p = mpatches.FancyArrow(0, 0.5*height, width, 0, length_includes_head=True, head_width=0.75*height )
    return p

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.plot([0,1,1,0,0],[0,0,0,0,0],[0,0,1,1,0],linewidth = 1.5, color = 'black')
plt.plot([1,2,2,1,1],[0,0,0,0,0],[0,0,1,1,0],linewidth = 1.5, color = 'black')
plt.plot([0,1,1,0,0],[0,0,0,0,0],[1,1,2,2,1],linewidth = 1.5, color = 'black')
plt.plot([1,2,2,1,1],[0,0,0,0,0],[1,1,2,2,1],linewidth = 1.5, color = 'black')
nodes = ax.scatter([0,1,2,0,1,2,0,1,2],[0,0,0,0,0,0,0,0,0],[0,0,0,1,1,1,2,2,2], s = 50, color = 'black', alpha = 1)
ax.quiver(1,0,1, 0.6667,0,0, pivot='tail', zorder = 4, color = 'C1')
ax.quiver(1,0,1, 0,0,0.6667, pivot='tail', zorder = 4, color = 'C1')
ax.quiver(1,0,1, 0,-0.6667,0, pivot='tail', zorder = 4, color = 'C0')


arr1 = plt.arrow(5,5,1,1, color = 'C0')
arr2 = plt.arrow(5,5,1,1, color = 'C1')

ax.set_xlim([0,2])
ax.set_ylim([0,2])
ax.set_zlim([0,2])
plt.axis('off')
ax.view_init(elev=12, azim=-47)
plt.legend([arr1,arr2,nodes],['Normal vector', 'Tangential vector', 'Sliding surface node'], handler_map={mpatches.FancyArrow : HandlerPatch(patch_func=make_legend_arrow),}, bbox_to_anchor=(0.6, .82), loc="upper left")
plt.savefig(figPath + "sliding_node_vecs_surf.png", dpi=250 ,bbox_inches='tight')
#%%
plt.close('all')
plt.rcParams["figure.figsize"] = [6.4, 4.8]
from matplotlib.legend_handler import HandlerPatch

def make_legend_arrow(legend, orig_handle,
                      xdescent, ydescent,
                      width, height, fontsize):
    p = mpatches.FancyArrow(0, 0.5*height, width, 0, length_includes_head=True, head_width=0.75*height )
    return p

pts = np.array([[-1.5,0,1.5],[0,0,0],[0,0,0]])

fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax = Axes3D(fig)
ax = fig.add_subplot(111, projection='3d')

ax.quiver(0,0,0, 0,0,1, pivot='tail', zorder = 4, label = '_nolegend_')
ax.quiver(0,0,0, 0,-1,0, pivot='tail', zorder = 4)
ax.quiver(0,0,0, 1,0,0, pivot='tail', zorder = 4, color = 'C1', label = 'Tangential vector')
arr1 = plt.arrow(5,5,1,1, color = 'C0')
arr2 = plt.arrow(5,5,1,1, color = 'C1')
ax.plot(pts[0,:], pts[1,:], pts[2,:], color = 'black', zorder = 3)
nodes = ax.scatter(pts[0,:], pts[1,:], pts[2,:], color = 'black', zorder = 5, s = 50, alpha = 1)
lim = 1.5
ax.set_xlim([-lim, lim])
ax.set_ylim([-lim, lim])
ax.set_zlim([-lim, lim])
plt.axis('off')
plt.legend([arr1,arr2,nodes],['Normal vector', 'Tangential vector', 'Sliding edge node'], handler_map={mpatches.FancyArrow : HandlerPatch(patch_func=make_legend_arrow),}, bbox_to_anchor=(0.6, .82), loc="upper left")
ax.view_init(elev=31, azim=-51)
plt.show()

plt.savefig(figPath + "sliding_node_vecs_3d.png", dpi=250 ,bbox_inches='tight')
#%%
plt.close('all')
fig = plt.figure()
ax = plt.gca()
nodes = plt.scatter([-.75,0,.75],[0,0,0],color = 'black',alpha = 1, s = 50, zorder = 5)
plt.plot([-.75,0,.75],[0,0,0],color = 'black', zorder = 3, linewidth = 1.5)
ax.quiver(0,0,0,-1, pivot = 'tail',scale= 6, color = 'C0', zorder= 4)
ax.quiver(0,0,1,0, pivot = 'tail',scale= 6, color = 'C1', zorder = 4)
arr1 = plt.arrow(5,5,1,1, color = 'C0')
arr2 = plt.arrow(5,5,1,1, color = 'C1')
plt.axis('equal')
plt.axis('off')
plt.xlim([-1.5,1.5])
plt.ylim([-1.4,1.6])
plt.legend([arr1,arr2,nodes],['Normal vector', 'Tangential vector', 'Sliding edge node'], handler_map={mpatches.FancyArrow : HandlerPatch(patch_func=make_legend_arrow),}, bbox_to_anchor=(0.6, .82), loc="upper left")
plt.show()
plt.savefig(figPath + "sliding_node_vecs_2d.png", dpi=250 ,bbox_inches='tight')


#%%
plt.close('all')
plt.rcParams["font.family"] = "Arial"    
plt.rcParams["font.size"] = 14
plt.rcParams["figure.figsize"] = [6.4, 0.8*4.8]
theta = np.linspace(0,np.pi,100)
x_circle = 1*np.cos(theta)+2.5
y_circle = 1*np.sin(theta)

x = np.array([0,5,5,0,0])
y = np.array([0,0,5,5,0])


# NON periodic
fig, ax = plt.subplots()
plt.scatter(2.5,0, s = 50, zorder = 2,label = 'Control node')
plt.plot(x,y,color='black', linewidth = 1.5, label = 'Domain boundary', zorder = 1)

plt.fill_between(x_circle,y_circle,y2=0, color = 'red',alpha = 0.5, label = "Compact RBF support")
plt.axis('equal')

plt.legend( bbox_to_anchor=(1.01, 1), loc="upper left")
plt.subplots_adjust(right=0.6)
plt.axis('off')
plt.show()
plt.savefig(figPath + "trans_per_support1.png", dpi=150 ,bbox_inches='tight')

# periodic
fig, ax = plt.subplots()
plt.scatter([2.5],[0], s = 50, zorder = 2, label = "Control node")
plt.plot(x,y,color='black', linewidth = 1.5, label = 'Domain boundary', zorder = 1)

yc = 0
xc = 2.5
l = 5
xr = np.linspace(1.5,3.5,100)
yr = 0 + l/np.pi * np.arcsin( np.pi/l * np.sqrt(1-(xr-xc)**2)   )
plt.fill_between(xr,yr,y2=0, color = 'red',alpha = 0.5, label = "Compact RBF support")

yr2 = 5-yr
plt.fill_between(xr,yr2,y2=5, color = 'red',alpha = 0.5, label = '_nolegend_')

plt.axis('equal')

plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
plt.subplots_adjust(right=0.6)
plt.axis('off')
plt.show()
plt.savefig(figPath + "trans_per_support2.png", dpi=150 ,bbox_inches='tight')
#%%
plt.rcParams["figure.figsize"] = [6.4, 0.4*4.8]
#Rotional periodic domain
r_outer = 10
theta = np.linspace(0,np.pi,100)
x_circle = 1*np.cos(theta)
y_circle = 1*np.sin(theta)

theta = np.linspace(0,np.pi/6)

theta2 = np.linspace(-np.pi*5/6,np.pi/6)
x_circle2 = 1*np.cos(theta2)    
y_circle2 = 1*np.sin(theta2)

x = np.hstack([0,r_outer*np.cos(theta),0])
y = np.hstack([0,r_outer*np.sin(theta),0])

plt.figure()
plt.plot(x,y,linewidth = 1.5, zorder = 1, color = 'black',label = 'Domain boundary')
plt.scatter(6.0,0, s = 50, zorder = 2, label = 'Control node')
plt.fill_between(x_circle+6.0,y_circle,y2=0, color = 'red',alpha = 0.5, label = "Compact RBF support")
plt.axis('equal')
plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
plt.subplots_adjust(right=0.6)
plt.axis('off')
plt.show()
plt.savefig(figPath + "rot_per_support1.png", dpi=150 ,bbox_inches='tight')
#%%

r_outer = 10
theta = np.linspace(0,np.pi,100)
x_circle = 1*np.cos(theta)
y_circle = 1*np.sin(theta)

theta = np.linspace(0,np.pi/6)

theta2 = np.linspace(-np.pi*5/6,np.pi/6)
x_circle2 = 1*np.cos(theta2)    
y_circle2 = 1*np.sin(theta2)

x = np.hstack([0,r_outer*np.cos(theta),0])
y = np.hstack([0,r_outer*np.sin(theta),0])

r_c = 6
t_c = 0

t_p = np.pi/6

r = np.linspace(5,7,100)
#print(np.pi/t_p * np.arccos( (-1+r**2+r_c**2)/(2*r*r_c)))

t = t_c + t_p/np.pi* np.arcsin( np.pi/t_p * np.arccos( (-(.85**2)+r**2+r_c**2)/(2*r*r_c) )  )
t2 = np.pi/6-t

x_s = r*np.cos(t)
y_s = r*np.sin(t)
x_s2 = r*np.cos(t2)
y_s2 = r*np.sin(t2)

x_s[:8] = 5.15
x_s[-8:] = 6.85
y_s[:8] = 0
y_s[-8:] = 0

x_s2[:8] = 4.46003082948986
x_s2[-8:] = 5.932274015923404
y_s2[:8] = 2.5749999999999997
y_s2[-8:] = 3.4249999999999994


plt.figure()
plt.plot(x,y,linewidth = 1.5, zorder = 1, color = 'black', label = "Domain boundary")
plt.scatter([6.0],[0], s = 50, zorder = 2, label = 'Control node')
plt.fill_between(x_s,y_s,y2=0, color = 'red',alpha = 0.5, label = "Compact RBF support")
y2 = 0.5773502690459182*x_s2
plt.fill_between(x_s2,y_s2,y2, color = 'red',alpha = 0.5, label = "_nolegend_")



plt.axis('equal')
plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left")
plt.subplots_adjust(right=0.6)
plt.axis('off')
plt.show()
plt.savefig(figPath + "rot_per_support2.png", dpi=150 ,bbox_inches='tight')
#%%
if(0):
    #%%
    plt.close('all')
    plt.rcParams["font.size"] = 16
    x = np.arange(0,5,0.01)
    
    phi = (1-abs(x))**4 * (4*abs(x)+1)
    phi_r2 = (1-abs(x/2))**4 * (4*abs(x/2)+1)
    phi_r5 = (1-abs(x/5))**4 * (4*abs(x/5)+1)
    phi[abs(x)>1] = 0
    phi_r2[abs(x)>2] = 0
    phi_r5[abs(x)>5] = 0
    
    plt.figure()
    plt.plot(x,phi,linewidth=2)
    plt.plot(x,phi_r2,linewidth=2)
    plt.plot(x,phi_r5,linewidth=2)
    plt.xlabel('$\\delta$ ')
    plt.ylabel('$\phi(\\xi)$')
    plt.legend(['$r=1$','$r=2$','$r=5$'])
    plt.savefig(figPath + "rbf_xi.png", dpi=800,bbox_inches = 'tight')
    #%%
    from matplotlib import cm
    from colMap import colMap
    from mpl_toolkits.mplot3d import Axes3D
    plt.close('all')
    plt.rcParams["font.family"] = "Arial"    
    plt.rcParams["figure.figsize"] = [6.4*1.3, 1.3*4.8]
    plt.rcParams["font.size"] = 16
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
    plt.xlabel('$x$')
    plt.ylabel('$y$')
    ax.zaxis.labelpad=20
    ax.yaxis.labelpad=10
    ax.xaxis.labelpad=10
    ax.zaxis.set_rotate_label(False)
    ax.set_zlabel('$\phi(\\xi)$', rotation = 0)
    v = np.linspace(0, 1.0, 6, endpoint=True)
    fig.colorbar(surf, shrink=0.5, ticks = v)
    ax.view_init(elev=37., azim=65)
    plt.savefig(figPath + "rbf_xy-plane.png", dpi=200)
    
    
    #%%
    plt.rcParams["font.size"] = 18
    x = np.linspace(0,1,100)
    l = 1
    y = l/np.pi*np.sin(x*np.pi/l)
    plt.figure()
    plt.plot(x,y, linewidth=2)
    
    
    plt.axis([0,1,0,1])
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlabel('$y-y_c$')
    plt.ylabel('$\delta_y$')
    plt.show()
    
    plt.savefig(figPath + "periodic_sine.png", dpi=150,bbox_inches='tight')
    
    
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
    g = plt.scatter(x[3],y[3],color = 'C0',zorder = 3)
    for i in range(len(n_x)-1):
    #i = 
        if i == 2:
            arrow1 = plt.arrow(midpntx[i], midpnty[i],n_x[i],n_y[i], head_width = .02, color='C1', zorder = 3, label = 'Edge normal vector')
    #    else:
    #        plt.arrow(midpntx[i], midpnty[i],n_x[i],n_y[i], head_width = .02, color='C1', zorder = 3, label = '_nolegend_')
    #    plt.arrow([midpntx[i], midpntx[i]+n_x[i]],[midpnty[i], midpnty[i]+n_y[i]], color='blue')
    #    plt.plot([x[i], x[i]+t_x[i]],[y[i], y[i]+t_y[i]], color='blue')
    d = plt.scatter(dx[2],dy[2], color = 'C2', zorder=3, label = 'Displaced sliding node')    
    
    xn = dx+5*n_x[:-1]*(deltax*(n_x[:-1]/np.sqrt(n_x[:-1]**2+n_y[:-1]**2)) + deltay*(n_y[:-1]/np.sqrt(n_x[:-1]**2+n_y[:-1]**2)))
    yn = dy+5*n_y[:-1]*(deltax*(n_x[:-1]/np.sqrt(n_x[:-1]**2+n_y[:-1]**2)) + deltay*(n_y[:-1]/np.sqrt(n_x[:-1]**2+n_y[:-1]**2)))
    #for i in range(len(deltax)):
    i =2
    e = plt.arrow(dx[i],dy[i], 3*n_x[i]*(deltax[i]*(n_x[i]/np.sqrt(n_x[i]**2+n_y[i]**2)) + deltay[i]*(n_y[i]/np.sqrt(n_x[i]**2+n_y[i]**2))), 3*n_y[i]*(deltax[i]*(n_x[i]/np.sqrt(n_x[i]**2+n_y[i]**2)) + deltay[i]*(n_y[i]/np.sqrt(n_x[i]**2+n_y[i]**2))),  head_width = .02, color = 'C2', label = 'Projection')
    #plt.scatter(xn,yn, color = 'green', zorder = 3)
    f = plt.scatter(xn[2],yn[2],color = 'C2', facecolor = 'none', zorder = 4)
    plt.axis('equal')
    plt.ylim((-0.1,1.25))
    plt.axis('off')
    
    plt.legend([b,c,arrow1,g,d,e, f], ['Sliding edge nodes','Nearest edge midpoint','Edge normal vector','Non-displaced sliding node','Displaced sliding node', 'Projection', 'Projected node'], handler_map={mpatches.FancyArrow : HandlerPatch(patch_func=make_legend_arrow),}, loc = 'lower right')
    plt.show()
    plt.savefig(figPath + "projection.png", dpi=200,bbox_inches='tight')
    
    
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
    g = plt.scatter(x[1],y[1], color = 'C0',zorder = 3)
    
    for i in range(len(n_x)):
    #i = 
        if i == 0:
            arrow1 = plt.arrow(mpx[i], mpy[i],-n_x[i],-n_y[i], head_width = .07, color='C1', zorder = 3, label = 'Line normal vector')
    #    else:
    #        plt.arrow(mpx[i], mpy[i],-n_x[i],-n_y[i], head_width = .07, color='C1', zorder = 3, label = '_nolegend_')
            
#    d = plt.scatter(-1.25,1.3,color = 'C2')
#    e = plt.arrow(-1.25,1.3,1.1,0, head_width = .07, color='C2', zorder = 3)
#    f = plt.scatter(0,1.3,color = 'C2', marker = 'o',facecolor='none')
    
    d = plt.scatter(0,1.3,color = 'C2')
    e = plt.arrow(0,1.3,0.05,-0.05, head_width = .07, color='C2', zorder = 3)
    f = plt.scatter(0.15,1.3-0.15,color = 'C2', marker = 'o',facecolor='none')
    plt.legend([b,c,arrow1,g,d,e,f], ['Sliding edge nodes','Nearest edge midpoint','Edge normal vector','Non-displaced sliding node','Displaced sliding node', 'Projection', 'Projected node'], handler_map={mpatches.FancyArrow : HandlerPatch(patch_func=make_legend_arrow),}, bbox_to_anchor=(0.5,0.71),loc = 'upper left')
    plt.axis('equal')
    plt.ylim((-1,3.))
    plt.xlim((-1.5,3))
    plt.axis('off')
    plt.show()
    plt.savefig(figPath + "projection3.png", dpi=200,bbox_inches='tight')