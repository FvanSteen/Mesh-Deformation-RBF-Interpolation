from matplotlib import pyplot as plt
from funs import getMeshQual, getPlotData
import numpy as np
from coordTransform import coordTransform
import os 

def stepPlots(f, s,p,g,de, save):
    
    figPath = os.path.dirname(os.path.abspath(__file__)) + "/figs/"
    max_steps = 0
    size_factor = 1.25
    plt.rcParams["font.family"] = "Arial"    
    plt.rcParams["font.size"] = 20
    plt.rcParams["axes.labelsize"] = 20
    plt.rcParams["figure.figsize"] = [6.4, 4.8*size_factor]
    plt.rcParams["legend.fontsize"] = 16
    
    start_fig = plt.gcf().number
    lw = 2.5
    ms = 10
    linestyle = ['-.','-','--',':','-.','-','--',':', '-.']
    markerstyle = ['o','^','d','s','8','p','h','D','p']
    c = ['C0','C1','C2','C3', 'C4','C5','C6','C7','C8']
    cnt = 0
    for fName in f:
        data = np.genfromtxt('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\runHistory\\'+fName+ ".txt", skip_header=1, dtype=None, encoding='utf-8')
        for smode in s:
            for pmode in p:
                for greedy in g:
                    for doubleEdge in de:
                        rows = ((data['f0']==smode) & (data['f1'] == pmode) & (data['f2'] == greedy) & (data['f3'] == doubleEdge)).nonzero()[0]
                        selected_data = data[rows]
                        steps = np.unique(selected_data['f6'])
                        
                        cpu_time = np.empty(np.size(steps),dtype = float)
                        q_min = np.empty(np.size(steps),dtype = float)
                        q_mean = np.empty(np.size(steps),dtype = float)
                        q_max = np.empty(np.size(steps),dtype = float)
                        
                        for i in range(len(steps)):
                            idx = (selected_data['f6'] == steps[i]).nonzero()[0]
                            cpu_time[i] = np.mean(selected_data['f7'][idx])/1000
                            
                            q_min[i] = selected_data['f8'][idx[0]]
                            q_mean[i] = selected_data['f9'][idx[0]]
                            q_max[i]= selected_data['f10'][idx[0]]
                            if q_min[i] < 0:
                                q_min[i] = np.nan
                                q_mean[i] = np.nan
                                q_max[i] = np.nan
                        
                        if(np.size(steps) != 0 and max_steps < np.max(steps)):
                            max_steps = np.max(steps)
#                        if cnt != 1:
                        plt.figure(start_fig)
                        plt.plot(steps,cpu_time,'-x', linestyle = linestyle[cnt], marker = markerstyle[cnt],markersize = ms, markerfacecolor = 'none', linewidth = lw, color = c[cnt])
                        
                        plt.figure(start_fig+1)
                        plt.plot(steps,q_min,'-x',linestyle = linestyle[cnt], marker = markerstyle[cnt],markersize = ms, markerfacecolor = 'none', linewidth = lw, color = c[cnt])
                        
                        plt.figure(start_fig+2)
                        plt.plot(steps,q_mean,'-x', linestyle = linestyle[cnt], marker = markerstyle[cnt],markersize = ms, markerfacecolor = 'none', linewidth = lw,color = c[cnt])
                    
#                        plt.figure(start_fig+3)
#                        plt.plot(cpu_time,q_min,'-x', linestyle = linestyle[cnt], marker = markerstyle[cnt],markersize = ms, markerfacecolor = 'none', linewidth = lw,color = c[cnt])
                        cnt += 1
                        
    plt.figure(start_fig)
    plt.xlabel('Number of steps [-]')
    plt.ylabel('CPU time [s]')
    plt.xticks(range(0,max_steps+1,2))
#    plt.legend(["Non-Per.","Per.","Per. Disp. (fixed)", "Per. Disp. (non-fixed)"],loc = 'upper left')
#    plt.legend(["Non-Per.","Per.","Per. Disp. (fixed)"],loc = 'upper left')
    plt.legend(["RBF", "RBF-PS","RBF-DS"],loc = 'lower right')
#    plt.legend(["RBF (fixed)","RBF-PS (fixed)","RBF-PS (non-fixed)","RBF-DS (fixed)","RBF-DS (non-fixed)"],bbox_to_anchor = [1,1],loc='upper left')
#    plt.legend(["RBF","RBF Per. Disp. (fixed)","RBF-PS","RBF-PS Per. Disp. (fixed)","RBF-DS","RBF-DS Per. Disp.  (fixed)"],bbox_to_anchor = [1,1],loc='upper left')
    plt.tight_layout()
    if save:
        plt.savefig(figPath + "25x25_moderate_def_cpu.png", dpi=100)
    
    plt.figure(start_fig+1)
    plt.xlabel('Number of steps [-]')
    plt.ylabel('Minimum qualitiy [-]')
    plt.xticks(range(0,max_steps+1,2))
#    plt.legend(["Non-Per.","Per.","Per. Disp. (fixed)", "Per. Disp. (non-fixed)"],loc='lower left')#,bbox_to_anchor=[0.35,.28])
#    plt.legend(["Non-Per.","Per.","Per. Disp. (fixed)"],loc = 'lower left')
    plt.legend(["RBF","RBF-PS","RBF-DS"],loc = 'lower right')#, bbox_to_anchor = [0.6,0.4])
#    plt.legend(["RBF (fixed)","RBF-PS (fixed)","RBF-PS (non-fixed)","RBF-DS (fixed)","RBF-DS (non-fixed)"],bbox_to_anchor = [1,1],loc='upper left')
#    plt.legend(["RBF","RBF Per. Disp. (fixed)","RBF-PS","RBF-PS Per. Disp. (fixed)","RBF-DS","RBF-DS Per. Disp.  (fixed)"],bbox_to_anchor = [1,1],loc='upper left')
#    plt.gca().set_ylim(bottom=0.2)
    plt.tight_layout()
    
    
    if save:
        plt.savefig(figPath + "25x25_moderate_def_qmin.png", dpi=100)
        
    plt.figure(start_fig+2)
    plt.xlabel('Number of steps [-]')
    plt.ylabel('Mean qualitiy [-]')
    plt.xticks(range(0,max_steps+1,2))
#    plt.legend(["Non-Per.","Per.","Per. Disp. (fixed)", "Per. Disp. (non-fixed)"],loc='lower left', bbox_to_anchor = [0.35,0.55])
#    plt.legend(["Non-Per.","Per.","Per. Disp. (fixed)"],loc = 'center right')
#    plt.legend(["RBF (fixed)","RBF-PS (fixed)","RBF-PS (non-fixed)","RBF-DS (fixed)","RBF-DS (non-fixed)"],bbox_to_anchor = [1,1],loc='upper left')
#    plt.legend(["RBF","RBF Per. Disp. (fixed)","RBF-PS","RBF-PS Per. Disp. (fixed)","RBF-DS","RBF-DS Per. Disp.  (fixed)"],bbox_to_anchor = [1,1],loc='upper left')
    plt.legend(["RBF","RBF-PS","RBF-DS"],loc = 'center right')
    
    plt.tight_layout()
    
    if save:
        plt.savefig(figPath + "25x25_moderate_def_qmean.png", dpi=100)
#        
#    plt.figure(start_fig+3)
#    plt.xlabel('CPU time [s]')
#    plt.ylabel('Minimum quality [-]')
##    plt.legend(["Non-Per.","Per.","Per. Disp. (fixed)"],loc = 'lower right')
##    plt.legend(["RBF (fixed)","RBF-PS (fixed)","RBF-PS (non-fixed)","RBF-DS (fixed)","RBF-DS (non-fixed)"],bbox_to_anchor = [1,1],loc='upper left')
#    plt.legend(["RBF","RBF-PS","RBF-DS"],loc = 'lower left')
#    plt.gca().set_ylim(bottom=0.2)
#    plt.tight_layout()
#    if save:
#        plt.savefig(figPath + "turbine_time_qmin.png", dpi=200)
    
    plt.show()

    
def defPlots(f, s,p,g,de, save, steps):
   
    figPath = os.path.dirname(os.path.abspath(__file__)) + "/figs/"

    size_factor = 1.25
    plt.rcParams["font.family"] = "Arial"    
    plt.rcParams["font.size"] = 20
    plt.rcParams["axes.labelsize"] = 20
    plt.rcParams["figure.figsize"] = [6.4, 4.8*size_factor]
    plt.rcParams["legend.fontsize"] = 16
#    start_fig = plt.gcf().number

    linestyle = ['-.','-','--',':','-.','-','--',':']
    markerstyle = ['o','^','d','s','o','^','d','s']
   
    q_min = np.empty((np.size(s),np.size(f)))
    q_mean = np.empty((np.size(s),np.size(f)))
    cpu_time = np.empty((np.size(s),np.size(f)))
    cnt = 0
    for fName in f:
        data =  np.genfromtxt('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\runHistory\\'+fName+ ".txt", skip_header=1, dtype=None, encoding='utf-8')
        cnt_s = 0
        for smode in s:
            
            for pmode in p:
                for greedy in g:
                    for doubleEdge in de:                        
                        rows = ((data['f0']==smode) & (data['f1'] == pmode) & (data['f2'] == greedy) & (data['f3'] == doubleEdge)).nonzero()[0]
                        selected_data = data[rows]    
                        idx = (selected_data['f6'] == steps).nonzero()[0]
                        
                        if(np.size(idx) != 0):         
                                                                           
                            q_min[cnt_s,cnt] = selected_data['f8'][idx[0]]
                            if q_min[cnt_s,cnt] < 0:
                                q_min[cnt_s,cnt] = np.nan
                            q_mean[cnt_s,cnt] = selected_data['f9'][idx[0]]
                            cpu_time[cnt_s,cnt] = np.mean(selected_data['f7'][idx])/1000
                            
                            
                        cnt_s += 1   
        cnt += 1
            
          
    
    plt.figure()
    for i in range(len(s)):
        plt.plot([0.0,0.2,0.4,0.6,0.8,1.0],cpu_time[i,:],linestyle = linestyle[i],marker = markerstyle[0], markerfacecolor = "none", linewidth = 1.5)

    plt.xlabel('Fraction of total deformation [-]')
    plt.ylabel('CPU time [s]')
    plt.legend(['RBF',"RBF-PS","RBF-DS"],loc = 'lower right')
    plt.gca().set_ylim(bottom=0)
    plt.tight_layout()
    if save:
        plt.savefig(figPath + "25x25_yVar_20_cpu.png", dpi=100)
        
        
    plt.figure()
    for i in range(len(s)):
        plt.plot([0.0,0.2,0.4,0.6,0.8,1.0],q_min[i,:],linestyle = linestyle[i],marker = markerstyle[0], markerfacecolor = "none", linewidth = 1.5)

    plt.xlabel('Fraction of total deformation [-]')
    plt.ylabel('Minimum quality [-]')
    plt.legend(['RBF',"RBF-PS","RBF-DS"],loc = 'upper right')
    plt.gca().set_ylim(bottom=0)
    plt.tight_layout()
    if save:
        plt.savefig(figPath + "25x25_yVar_20_qmin.png", dpi=100)
        
    
    plt.figure()
    for i in range(len(s)):
        plt.plot([0.0,0.2,0.4,0.6,0.8,1.0],q_mean[i,:],linestyle = linestyle[i],marker = markerstyle[0], markerfacecolor = "none", linewidth = 1.5)

    plt.xlabel('Fraction of total deformation [-]')
    plt.ylabel('Mean quality [-]')
    plt.legend(['RBF',"RBF-PS","RBF-DS"],loc = 'upper right')
    plt.gca().set_ylim(bottom=0)
    plt.tight_layout()
    if save:
        plt.savefig(figPath + "25x25_yVar_20_qmean.png", dpi=100)
        
#    plt.figure()
#    for i in range(len(s)):
#        plt.plot(cpu_time[i,:],q_min[i,:],linestyle = linestyle[i],marker = markerstyle[0], markerfacecolor = "none", linewidth = 1.5)
#
#    plt.xlabel('Fraction of total deformation')
#    plt.ylabel('Mean quality [-]')
#    plt.legend(['RBF',"RBF-PS","RBF-DS"],loc = 'upper right')
#    plt.gca().set_ylim(bottom=0)
#
#    if save:
#        plt.savefig(figPath + fName + "_time_qmin.png", dpi=400)


def qualDistribution(fileNames, cutAxis, labels, qType):
    plt.rcParams["font.family"] = "Arial"    
    plt.rcParams["font.size"] = 14
    dims = [0,1,2]
    dims.remove(cutAxis)
    
    
    for i in range(np.size(fileNames)):
        [f,v,_,_,_,_,_] = getPlotData(fileNames[i]) 
        v = coordTransform.toCylindrical(v)
        meshQuality = getMeshQual(fileNames[i])

        maxVal = np.max(v[:,cutAxis])
        minVal = np.min(v[:,cutAxis])
        
        cutPlanes = np.linspace(minVal,maxVal,102)

        cutPlanes = np.delete(cutPlanes,0)
        cutPlanes = np.delete(cutPlanes,-1)
        
        if(i == 0):
            minQ = np.empty((np.size(fileNames),np.size(cutPlanes)))
        
        for j in range(np.size(cutPlanes)):
            
            idx_lower = np.where(v[f][:,:,cutAxis] < cutPlanes[j])
            idx_upper = np.where(v[f][:,:,cutAxis] >= cutPlanes[j])
            idxCutElems = np.intersect1d(idx_lower[0],idx_upper[0])
            
            if(np.size(idxCutElems)==0):
                minQ[i,j] = np.nan
            else:
                meshQual_cut = meshQuality[idxCutElems]
                if qType == "min":
                    minQ[i,j] = np.min(meshQual_cut)
                elif qType == "mean":
                    minQ[i,j] = np.mean(meshQual_cut)
                else:
                    print("Invaled quality type!")
                    quit
                    

    plt.figure()
    for row in range(np.size(fileNames)):    
        plt.plot(cutPlanes,minQ[row,:], lineWidth = 2)
   
    if qType == "min":
        plt.ylabel("Min. mesh quality [-]")
    elif qType == "mean":
        plt.ylabel("Mean. mesh quality [-]")
    
    
    if cutAxis == 0:
        plt.xlabel("r")
    elif cutAxis == 1:
        plt.xlabel("$\\theta$")
    elif cutAxis == 2:
        plt.xlabel("z")
    plt.legend(labels,fontsize = 12)
#    plt.ylim([0,1.02])
    plt.show()
    
        