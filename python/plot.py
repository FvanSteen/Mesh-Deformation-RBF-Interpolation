import os
from bdryScatter import bdryScatterFuns
from meshQualPlot import meshQualPlot2D, meshQualPlot3D, nanElems
from funs import getMeshQualParams,getPlotData, getMeshQualParams3D, getMeshQuals3D, getMeshQual, getMeshQuals
import matplotlib.pyplot as plt
import numpy as np
from coordTransform import coordTransform
from ThesisDataPlots import stepPlots, defPlots, qualDistribution

import matplotlib.collections
#plt.close('all')
# Setting directory to find .su2 files
figPath = os.path.dirname(os.path.abspath(__file__)) + "/figs/"
os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\Meshes')

#%%


#fig = plt.figure()
#ax = fig.add_subplot(1,1,1, projection='3d')
#
#elem = v[f][116049]
#elem2 = elem
#print(elem)
#verts = [np.array([elem[0],elem[1],elem[5],elem[4]]),
#      np.array([elem[3],elem[2],elem[6],elem[7]]),
#      np.array([elem[3],elem[0],elem[4],elem[7]]),
#      np.array([elem[2],elem[1],elem[5],elem[6]]),
#      np.array([elem[0],elem[1],elem[2],elem[3]]),
#      np.array([elem[4],elem[5],elem[6],elem[7]])]
#pc = Poly3DCollection(verts, facecolors='red', edgecolor="black",linewidth=0.5)
#    
#ax.add_collection(pc)
#
#ub = np.max(elem,axis=0)
#lb = np.min(elem,axis=0)
#ax.set_xlim(lb[0], ub[0])
#ax.set_ylim(lb[1], ub[1])
#ax.set_zlim(lb[2], ub[2])
#
#fileName = '/stator_per_ffd_mod_def.su2'
##fileName = '/25x25_per_def.su2'
#[f,v,elemType,bdryPnts,markerTags, nElemsMarks, FFD_pnts] = getPlotData(fileName)  
#elem = v[f][116049]
#elem1 = elem
#print(elem)
#verts = [np.array([elem[0],elem[1],elem[5],elem[4]]),
#      np.array([elem[3],elem[2],elem[6],elem[7]]),
#      np.array([elem[3],elem[0],elem[4],elem[7]]),
#      np.array([elem[2],elem[1],elem[5],elem[6]]),
#      np.array([elem[0],elem[1],elem[2],elem[3]]),
#      np.array([elem[4],elem[5],elem[6],elem[7]])]
#pc = Poly3DCollection(verts, facecolors='none', edgecolor="black",linewidth=0.5)
#    
#ax.add_collection(pc)
#
#ub = np.max(elem,axis=0)
#lb = np.min(elem,axis=0)
#ax.set_xlim(lb[0], ub[0])
#ax.set_ylim(lb[1], ub[1])
#ax.set_zlim(lb[2], ub[2])

#
#plt.close('all')
#idx = [322,327,353,348,322]
#l = np.arange(0,26,1)
#ll = np.arange(25,676,26)
#u = np.arange(650,676,1)
#r = np.arange(0,651,26)
#plt.figure()
##plt.plot([0,1,1,0,0],[0,0,1,1,0],linewidth = 1.0,color = 'black', linestyle = '-',label='_nolegend_')
#plt.plot([0.4,0.6,0.6,0.4,0.4],[0.48,0.48,0.52,0.52,0.48],linewidth = 1.5, linestyle = '-', color = 'black')
#plt.plot(v[l][:,0],v[l][:,1],linewidth = 1.0,color = 'black', linestyle = '-',label='_nolegend_')
#plt.plot(v[ll][:,0],v[ll][:,1],linewidth = 1.0,color = 'black', linestyle = '-',label='_nolegend_')
#plt.plot(v[u][:,0],v[u][:,1],linewidth = 1.0,color = 'black', linestyle = '-',label='_nolegend_')
#plt.plot(v[r][:,0],v[r][:,1],linewidth = 1.0,color = 'black', linestyle = '-',label='_nolegend_')
#plt.plot(v[idx][:,0],v[idx][:,1],linewidth = 1.5, linestyle = '--', color = 'black')
##fileName = '/25x25_per_def.su2'
##[f,v,elemType,bdryPnts,markerTags, nElemsMarks, FFD_pnts] = getPlotData(fileName) 
#plt.plot(v[idx][:,0],v[idx][:,1],linewidth = 1.5, linestyle = '--', color = 'black')
#plt.legend(['Initial\nposition','Deformed\nposition'], bbox_to_anchor=(1.0,1), loc='upper left')
#
#plt.tight_layout()
#plt.axis([-0.05, 1.05, -0.05, 1.05])
#plt.gca().set_aspect('equal', adjustable='box')
#plt.show()


#plt.close('all')
#
#init_x = np.array([0.4,0.6,0.6,0.4,0.4])
#init_y = np.array([0.48,0.48,0.52,0.52,0.48])
#
#dx = -0.35/5
#dy = -0.43/5
#
#plt.figure()
#plt.plot([0,1,1,0,0],[0,0,1,1,0],linewidth = 1.0,color = 'black', linestyle = '-',label='_nolegend_')
#plt.plot(init_x,init_y,linewidth = 1.5, linestyle = '-', color = 'black')
#plt.plot(init_x+dx,init_y+dy,linewidth = 1.5, linestyle = '--', color = 'black')
#plt.text(init_x[1]+dx+0.025,init_y[1]+1*dy+0.005, "Loc 1")
#plt.plot(init_x+2*dx,init_y+2*dy,linewidth = 1.5, linestyle = '--', color = 'black')
#plt.text(init_x[1]+2*dx+0.025,init_y[1]+2*dy+0.005, "Loc 2")
#plt.plot(init_x+3*dx,init_y+3*dy,linewidth = 1.5, linestyle = '--', color = 'black')
#plt.text(init_x[1]+3*dx+0.025,init_y[1]+3*dy+0.005, "Loc 3")
#plt.plot(init_x+4*dx,init_y+4*dy,linewidth = 1.5, linestyle = '--', color = 'black')
#plt.text(init_x[1]+4*dx+0.025,init_y[1]+4*dy+0.005, "Loc 4")
#plt.plot(init_x+5*dx,init_y+5*dy,linewidth = 1.5, linestyle = '--', color = 'black')
#plt.text(init_x[1]+5*dx+0.025,init_y[1]+5*dy+0.005, "Loc 5")
#plt.legend(['Initial\nposition','Deformed\nposition'], bbox_to_anchor=(1.04,1), loc='upper left')
#plt.tight_layout()
#plt.axis([-0.05, 1.05, -0.05, 1.05])
#plt.gca().set_aspect('equal', adjustable='box')
#
#plt.show()
##
#plt.savefig(figPath + "25x25_square_deformation.png", dpi=300,bbox_inches='tight')

#v = coordTransform.toPolar(v)
#%% 2D PLOTS
#plt.close("all")
# Mesh quality

#from mpl_toolkits.mplot3d.art3d import Poly3DCollection
#fileName = '/25x25x25_per_def.su2'
#fileName = '/stator_per_ffd_ps.su2'
#fileName = '/nthx001_identical_tube4_mesh_def.SU2'
#fileName = '/200x200_def.su2'
#fileName = '/turbine_row.su2'
for i in range(20):
    fileName = '/25x25_def'+str(i)+'.su2'
#    fileName = '/25x25_def.su2'
    [f,v,elemType,bdryPnts,markerTags, nElemsMarks, FFD_pnts] = getPlotData(fileName)
#    v = coordTransform.toCylindrical(v)
    
    
    plt.rcParams["font.family"] = "Arial"
    plt.rcParams["font.size"] = 16
    
    if(1):    
        graphName = ''   
        plt.rcParams["font.size"] = 16
        meshQualPlot2D(fileName, graphName, f, v, elemType)
        plt.savefig(figPath + "sq_ds_"+str(i)+".png", dpi=200,bbox_inches='tight')
#        plt.savefig(figPath + "sq_moving.png", dpi=800,bbox_inches='tight')
        plt.close("all")
# Boundary scatter plots
if(0):
    # Without initial boundary
    init = 0
    plt.rcParams["font.size"] = 16
    plt.rcParams["figure.figsize"] = [6.4, 4.8]
    bdryScatterFuns.bdryScatter(fileName, v, bdryPnts, init)
    
#    plt.legend(['Initial\n position', 'Deformed\nposition', 'Enlarging\nreference'], bbox_to_anchor = [1,0.6], loc = 'lower left')
#    plt.savefig(figPath + "turbine_deformed.png", dpi=200,bbox_inches='tight')
    plt.savefig(figPath + "per3.png", dpi=800, bbox_inches='tight')
    

    

        
#%%

#ax.legend(labels, bbox_to_anchor=(1.04,1), loc='upper left')
#plt.show()
#plt.savefig(figPath + "rbf_def.png", dpi=800,bbox_inches='tight')


#%% 3D PLOTS
    
if(0):
    cutAxis = 0
#    cutPlaneLoc = 0.2
    
    fileNames = ['/stator_per_ffd_base.su2','/stator_per_ffd_elastic_deform.su2', '/stator_per_ffd_ps.su2','/stator_per_ffd_ds.su2']#,'/stator_per_ffd_test_def.su2']
#    fileNames = ['/stator_per_ffd_test_def.su2']
    labels = ['Initial','Elastic eqs.','RBF-PS','RBF-DS']
    qType = "min"
    qualDistribution(fileNames, cutAxis,labels, qType)
    plt.savefig(figPath + "aachen_radial_min_comp.png", dpi=400,bbox_inches='tight')
    
    
# Mesh quality
if(0):
#    plt.rcParams["figure.figsize"] = [6.4, 4.8]
    graphName = ''
    cutAxis = 0
    cutPlaneLoc = 0.2512
#    cutPlaneLoc = 0.6
    plt.rcParams["font.size"] = 14
#    plt.rcParams["figure.figsize"] = [6.4, 4.8]
    meshQualPlot3D(fileName, graphName, cutAxis, cutPlaneLoc, f, v, elemType)
    plt.savefig(figPath + "stator_r=0.2512_ps_double.png", dpi=300,bbox_inches='tight')
    
# Boundary scatter
if(0):
#    plt.rcParams["figure.figsize"] = [6.4, 4.8/1.5]
    plt.rcParams["font.size"] = 16
    plotTag = ["BLADE","HUB", "SHROUD","INFLOW","OUTFLOW"]
#    plotTag = ["HUB","PER2"]
#    plotTag = ["BLOCK","LEFT","RIGHT","UPPER","LOWER","FRONT","BACK"]
#    plotTag = ["LEFT", "RIGHT", "UPPER","LOWER","FRONT","BACK","BLOCK"]
    bdryScatterFuns.bdryScatter3D(v, bdryPnts, markerTags, nElemsMarks, plotTag, FFD_pnts)
    plt.savefig(figPath + "stator_node_div_hd_loose.png", dpi=800)

if(0):
    nanElems(fileName,v,f)
#    largerOneElems(fileName,v,f)
    
#%%
if(0):
    plt.close('all')
    fName = ["25x25_moderate_def"]
    s = ['ps']
    p = ['none', 'periodic', 'fixed','moving']
    g = [0] 
    de = [0]
    
    save_fig= False
    stepPlots(fName,s,p,g,de, save_fig)
    
if (0):
    fName = ["25x25_y0_def", "25x25_y1_def","25x25_y2_def","25x25_y3_def","25x25_y4_def","25x25_y5_def"]
    s = ['none','ps','ds']
    p = ['none']
    g = [0]
    de = [0]
    steps = 20
    save_fig= True
    defPlots(fName,s,p,g,de, save_fig, steps)
    
if(0):
    #fNames = [ "200x200.su2_ML_DE_size_128","200x200.su2_ML_DE_size_64","200x200.su2_ML_DE_size_32","200x200.su2_ML_DE_tol_0.1","200x200.su2_ML_DE_tol_0.2" ]#,"200x200.su2_SL_DE"]
    #fNames = ["200x200.su2_ML_DE_size_64", "200x200.su2_ML_DE_size_128","200x200.su2_ML_DE_size_256","200x200.su2_ML_DE_size_512","200x200.su2_SL_DE",]#, "200x200.su2_ML_DE_tol_0.1", "200x200.su2_ML_DE_tol_0.05" ]
    #fNames = ["200x200.su2_SL_SE","200x200.su2_SL_DE", "200x200.su2_ML_DE_tol_0.5", "200x200.su2_ML_DE_tol_0.2", "200x200.su2_ML_DE_tol_0.1", "200x200.su2_ML_DE_tol_0.05"]
    #fNames = ["200x200.su2_ML_DE_size_512","200x200.su2_ML_SE_size_512","200x200.su2_ML_DE_size_256","200x200.su2_ML_SE_size_256","200x200.su2_ML_DE_size_128","200x200.su2_ML_SE_size_128"]
    plt.rcParams["figure.figsize"] = [6.4, 4.8]
    plt.rcParams["axes.labelsize"] = 18
    plt.rcParams["font.size"] = 16
    plt.rcParams["legend.fontsize"] = 13
#    fNames = ["200x200.su2_ML_DE_size_8","200x200.su2_ML_DE_size_16","200x200.su2_ML_DE_size_32","200x200.su2_ML_DE_size_64", "200x200.su2_SL_de"]#,"200x200.su2_ML_DE_size_64","200x200.su2_ML_DE_size_128","200x200.su2_ML_DE_size_256", "200x200.su2_SL_DE"]#,"200x200.su2_ML_SE_size_128","200x200.su2_ML_DE_size_128", "200x200.su2_ML_SE_size_256","200x200.su2_ML_DE_size_256"]
#    fNames = ["stator_per_ffd_mod.su2_SL_DE_ps","stator_per_ffd_mod.su2_SL_DE_ds"]
    fNames = ["200x200.su2_ML_DE_tol_0.75","200x200.su2_ML_DE_tol_0.5","200x200.su2_ML_DE_tol_0.2","200x200.su2_ML_DE_tol_0.1","200x200.su2_ML_DE_tol_0.05","200x200.su2_SL_de"]#,"200x200.su2_ML_DE_tol_0.1", "200x200.su2_ML_DE_tol_0.05", "200x200.su2_SL_DE"]
#    fNames = ["200x200.su2_ML_DE_size_32","200x200.su2_ML_DE_size_64", "200x200.su2_ML_DE_size_128", "200x200.su2_ML_DE_size_256", "200x200.su2_SL_DE"]
    plt.close('all')
    for fName in fNames:
        data =  np.genfromtxt('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\convHist\\'+fName+ ".txt", skip_header=1, dtype=None, encoding='utf-8')
        t = data['f3']
        err = data['f2']
        N = data['f4']
        chord = 0.2#0.062128
        plt.figure(1)
        plt.plot(t/1000,err/chord, linewidth = 2, alpha=.65)
    #    plt.figure(2)
    #    plt.plot(t/1000, N, linewidth = 1.5)
    plt.figure(1)
    plt.yscale('log')
    plt.xlabel('CPU time [s]')
    plt.xlim([0,50])
    plt.ylabel(r'$\frac{\epsilon_{max}}{c}$ [-]')
    plt.legend(["S=8","S=16","S=32","S=64","Single level"])
#    plt.legend(["RBF-PS","RBF-DS"])
#    plt.legend(["Single edge","Double edge"])
#    plt.legend(["S=32","S=64","S=128","S=256","Single level"])
    plt.legend([r"$a_\epsilon$ = 0.75", r"$a_\epsilon$ = 0.5",r"$a_\epsilon$ = 0.2",r"$a_\epsilon$ = 0.1",r"$a_\epsilon$ = 0.05","Single level"])
    plt.savefig(figPath + "Multi_tol_compare_ds.png", dpi=150, bbox_inches = 'tight')
    plt.tight_layout()
    
    #plt.figure(2)
    ##plt.yscale('log')
    #plt.xlabel('CPU time [s]')
    #plt.xlim([0,20])
    #plt.ylabel('Number of control nodes')
    ##plt.legend(["128","256","512","SL-DE"])
    #plt.legend(["SE","DE"])
    #plt.tight_layout()
    
    plt.show()
 #%%
#from mpl_toolkits.mplot3d.art3d import Poly3DCollection
#fileName = '/stator_per_ffd_mod.su2'
##fileName = '/5x5x5_per_def.su2'
#[f_init,v_init,_,_,_,_,_] = getPlotData(fileName)  
#elem = v_init[f_init][76298]
#
#print(elem)
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1, projection='3d')
#
#verts = [np.array([elem[0],elem[1],elem[5],elem[4]]),
#      np.array([elem[3],elem[2],elem[6],elem[7]]),
#      np.array([elem[3],elem[0],elem[4],elem[7]]),
#      np.array([elem[2],elem[1],elem[5],elem[6]]),
#      np.array([elem[0],elem[1],elem[2],elem[3]]),
#      np.array([elem[4],elem[5],elem[6],elem[7]])]
#        
#pc = Poly3DCollection(verts, facecolors='red', edgecolor="black")
#ax.add_collection(pc)
#elem = v[f][76298]
#print(elem)
#verts = [np.array([elem[0],elem[1],elem[5],elem[4]]),
#      np.array([elem[3],elem[2],elem[6],elem[7]]),
#      np.array([elem[3],elem[0],elem[4],elem[7]]),
#      np.array([elem[2],elem[1],elem[5],elem[6]]),
#      np.array([elem[0],elem[1],elem[2],elem[3]]),
#      np.array([elem[4],elem[5],elem[6],elem[7]])]
#        
#pc = Poly3DCollection(verts, facecolors='blue', edgecolor="black")
#    
#ax.add_collection(pc)
#ax.set_xlim(0.290,0.303)
#ax.set_ylim(-0.065,-0.042)
#ax.set_zlim(0.0398,0.0408)
#plt.show()
#%%
if(0):
    #import numpy as np
    [f_init,v_init,elemType,bdryPnts_init,markerTags, nElemsMarks, FFD_pnts] = getPlotData("/stator_per_ffd_def.su2")
#    [f_init,v_init,elemType,bdryPnts_init,markerTags, nElemsMarks, FFD_pnts] = getPlotData("/OneraM6_fixed.su2")
#    plotTag= ["LOWER", "UPPER","LEFT","RIGHT","FRONT","BACK", "BLOCK"]
#    plotTag = ["BLADE", "HUB", "SHROUD"]
    plotTag = ["BLADE","HUB","INFLOW","OUTFLOW","PER1","PER2","SHROUD"]
#    plotTag = ["LOWER_SIDE" ,"UPPER_SIDE","TIP", "XNORMAL_FACES", "YNORMAL_FACES", "ZNORMAL_FACES", "SYMMETRY_FACE"]
    bdryScatterFuns.bdryScatter3D(v_init, bdryPnts_init, markerTags, nElemsMarks, plotTag, FFD_pnts)


#%%
#pnts = []
#for i in range(np.shape(v_init)[0]): 
#    if v_init[i,0] < 1.5 and v_init[i,0] > -0.5 and v_init[i,2] < 0.25 and v_init[i,2] > -0.25 and v_init[i,1] >= 0 and v_init[i,1] <1.3:
#        pnts.append(i)
#    
#fig = plt.figure()
#ax = fig.add_subplot(projection='3d')
#ax.scatter(v_init[pnts,0],v_init[pnts,1],v_init[pnts,2])
# 
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.set_zlabel('z')
#%%
#from mpl_toolkits import mplot3d
#
#fig = plt.figure()
#ax = plt.axes(projection="3d")

#ax.scatte3D()
#fname = '/25x25x5_def.su2'
#bdryScatterFuns.bdryScatter3D(fname)
#%%
if(0):
    import numpy as np
    import math
    import matplotlib.collections
    from colMap import colMap
    
    cmapMatlab = colMap()
#    fileNames = ['/OneraM6_fixed.su2']
    fileNames = ['/stator_per_ffd_def.su2']
#    fileNames = ['/5x5x5_per.su2']
    
#    [f_init,v_init,elemType,bdryPnts_init,_,_,_] = getPlotData('/OneraM6_fixed.su2')
    
    
    graphNames = ["",""]
    cutAxis = 0
    cutPlaneLoc = 0.27
    
    dims = [0,1,2]
    dims.remove(cutAxis)
    [f,v,elemType,_, markerTags, nElemsMarks,FFD_pnts] = getPlotData(fileNames[0]) 
       
    
    idx_lower = np.where(v[f][:,:,cutAxis] < cutPlaneLoc)
    idx_upper = np.where(v[f][:,:,cutAxis] >= cutPlaneLoc)
    idxCutElems = np.intersect1d(idx_lower[0],idx_upper[0])
    
    meshQuality = getMeshQual(fileNames[0])
    meshQual_cut = meshQuality[idxCutElems]
    

    print("Min mesh quality: \t", round(np.min(meshQual_cut),5)) 
    print("Max mesh quality: \t", round(np.max(meshQual_cut),5)) 
    print("Mean mesh quality: \t", round(np.mean(meshQual_cut),5))
    
    # quadriliteral
    verticesInd = np.array([[0,1],[1,2],[2,3],[3,0],[0,4],[1,5],[2,6],[3,7],[4,5],[5,6],[6,7],[7,4]])
    
    # tetrahidral
#    verticesInd = np.array([[0,1],[1,2],[2,0],[0,3],[1,3],[2,3]])
    
    v_cut = np.empty((len(idxCutElems),8,2),dtype = float)
    v_cut[:] = np.nan    
    
    for i in range(len(idxCutElems)):
        if(i%100 == 0):
            print(i)
#    for i in range(1):
        
        edges = v[f][idxCutElems[i]][verticesInd]   
        
#        cutLoc = np.empty((12,2), dtype = float)
#        cnt = 0
#        for ii in range(len(edges)):
#            if np.any(edges[ii][:,cutAxis] < cutPlaneLoc) and np.any(edges[ii][:,cutAxis] >= cutPlaneLoc):
#                delta = edges[ii][1,:] - edges[ii][0,:]
#                cutLoc[cnt,:] = edges[ii][0,dims] + delta[dims]*((cutPlaneLoc-edges[ii][0,cutAxis])/delta[cutAxis])
#                cnt = cnt+1
#                
        idx_1 = np.where(edges[:,:,cutAxis] < cutPlaneLoc)
        idx_2 = np.where(edges[:,:,cutAxis] >= cutPlaneLoc)
        idx_crossing = np.intersect1d(idx_1[0],idx_2[0])
        
        cutLoc = edges[idx_crossing][:,0,dims] + (edges[idx_crossing][:,1,dims]-edges[idx_crossing][:,0,dims])*((cutPlaneLoc-edges[idx_crossing][:,0,cutAxis])/(edges[idx_crossing][:,1,cutAxis]-edges[idx_crossing][:,0,cutAxis]))[:, np.newaxis]
        
        
        midpoint = np.mean(cutLoc,axis=0)
        
        
        
#        angles = np.empty(cnt,dtype=float)
#        for ii in range(cnt):
#            angles[ii] = math.atan2(cutLoc[ii,1]-midpoint[1],cutLoc[ii,0]-midpoint[0])
#    
#        print(angles)
        angles = np.arctan2(cutLoc[:,1]-midpoint[1],cutLoc[:,0]-midpoint[0])
        cutLoc = cutLoc[angles.argsort()]
        
        v_cut[i,0:np.shape(cutLoc)[0],:] = cutLoc
    #    print(cutLoc[0:cnt,:])
    
    #    plt.scatter(cutLoc[0:cnt,0],cutLoc[0:cnt,1], color= "blue") 
                
    #    for ii in range(len(v1)):
                
    #        print(v[f][idxCutElems[i]][v1,cutAxis], v[f][idxCutElems[i]][v2,cutAxis])
     

    # plotting stuff   
    colors = cmapMatlab(plt.Normalize(0,1)(meshQual_cut))
    fig = plt.figure()  
    pc = matplotlib.collections.PolyCollection(v_cut,cmap=cmapMatlab, facecolors=colors, edgecolor="black",linewidth=0.1)
    ax = fig.add_subplot(1,1,1)
    
    polys = ax.add_collection(pc)
    pc.set_array(None)
    ax.autoscale()  
    ax.set_aspect('equal')
    polys.set_clim(0,1)
    plt.colorbar(polys, ax=ax, shrink=1.0/(np.size(graphNames)-1))
    plt.show()
#%%
if(0):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    cutElems = v[f[idxCutElems]]
    for i in range(len(cutElems)):
        ax.scatter(cutElems[i][:,0],cutElems[i][:,1],cutElems[i][:,2],marker = "1", color='blue')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_zlim(0,1)
    #plotTag = ["BLOCK", "LOWER","UPPER","FRONT","BACK","LEFT","RIGHT"]
    #bdryScatterFuns.bdryScatter3D(v, bdryPnts_init, markerTags, nElemsMarks, plotTag, FFD_pnts)
    

#%%
#import os 
#import numpy as np
#os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\Classical RBF\\defs')
#x = os.getcwd()
#o = open(x + "\\gis_nthx001_mesh_deformation.txt", "w")
#i_base = open(x + "\\MoveSurface_baseline.txt", "r") 
#i_def = open(x + "\\MoveSurface_deformed.txt", "r") 
#
#lines_base = i_base.read().splitlines()
#lines_def = i_def.read().splitlines()
#
#delta = np.empty((len(lines_base), 2), dtype = float)
#index = np.empty((len(lines_base),1), dtype = int)
#for i in range(0,len(lines_base)):
#    i_split_base = lines_base[i].split()
#    i_split_def = lines_def[i].split()
#    delta[i,0] = float(i_split_def[1]) - float(i_split_base[1])
#    delta[i,1] = float(i_split_def[2]) - float(i_split_base[2])
#    index[i] = int(i_split_def[0])
##    print(i_split_base)
##    delta[i,0] = i_split[]
##    print(i.split()[0])
#
#
#for idx, c1,c2 in zip(index[:,0],delta[:,0],delta[:,1]):
#    o.write("{0}\t{1}\t{2}".format(idx, c1,c2) + "\n")
#o.close()

