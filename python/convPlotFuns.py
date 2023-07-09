import os
import numpy as np
import matplotlib.pyplot as plt

def convPlot(files, tol):
    plt.figure()
    xValMax = 0 
    for i in range(len(files)):
        startIdx = files[i].rfind("m")
        endIdx = files[i].rfind(".txt")
     
        [time, maxError,_] = getConvPlotData(files[i])
        xValMax = max(xValMax, max(time))
        plt.plot(time, maxError, label= files[i][startIdx+1:endIdx], marker = ".")
        
#    plt.plot([0,xValMax],[tol, tol], linestyle= '--')
        
    plt.xlabel('Computational time [s]')
    plt.ylabel('max error') # make non dimensional units
    plt.legend()
    plt.yscale('log')
    plt.show()
    
    
def getConvPlotData(fname):
    x = os.getcwd()
    
    fileObj = open(x+"\\convHist\\"+fname, "r") 
 
    lines = fileObj.read().splitlines()
    
    maxError = np.empty((len(lines)-1,1))    
    time = np.empty((len(lines)-1,1))    
    N = np.empty((len(lines)-1,1))    
    
    
    for i in range(1,len(lines)):
        data = lines[i].split()
        maxError[i-1] = float(data[2])
        time[i-1] = float(data[4])
        N[i-1] = int(data[5])
    return([time, maxError, N]) 
    
    
def plotCompTime(files):

    plt.figure()
    for file in files:
        [time,_,N] = getConvPlotData(file)
          
        
        dt = np.empty((len(time),1))
        dt[0] = time[0]
        for i in range(len(time)-1):
            dt[i+1] = time[i+1]-time[i]
        
        startIdx = file.rfind("m")
        endIdx = file.rfind(".txt")          
          
        if int(file[startIdx+1:endIdx]) != 1:
            for i in range(len(time)-1):
                N[i+1] = N[i+1]+N[i]
          
         
              
     
        
        plt.plot(N, dt,label= "M"+file[startIdx+1:endIdx], marker = ".")
    plt.legend()
    plt.ylabel('dt [s]')
    plt.xlabel('number control nodes') # make non dimensional units
#    
os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool')
#files = ["5x5_m1.txt", "5x5_m4.txt", "5x5_m8.txt","5x5_m12.txt", "5x5_m16.txt"]
#files = ["25x25.su2_ML_DE_tol0.1.txt","25x25.su2_SL_DE.txt","25x25.su2_ML_DE_tol0.2.txt"]
#files = ["375x375_m32.txt","375x375_m16.txt","375x375_m1.txt"]#,"375x375_m36.txt","375x375_m48.txt"]
#files = ["mesh_NACA0012_inv_m16.txt","mesh_NACA0012_inv_m8.txt","mesh_NACA0012_inv_m4.txt","mesh_NACA0012_inv_m1.txt"]#,"mesh_NACA0012_inv_m4.txt","mesh_NACA0012_inv_m1.txt"]
#files = ["LS89_turbine_m1.txt","LS89_turbine_m8.txt","LS89_turbine_m16.txt" ]
#files = ["OneraM6_fixed.su2_SL_DE.txt","OneraM6_fixed.su2_SL_SE.txt","OneraM6_fixed.su2_ML_DE_tol0.1.txt"]
#files = ["OneraM6_fixed.su2_ML_DE_tol0.1_r1.txt","OneraM6_fixed.su2_ML_DE_tol0.1_r2.txt","OneraM6_fixed.su2_ML_DE_tol0.1.txt"]
#files = ["OneraM6_fixed.su2_ML_DE_tol0.5.txt","OneraM6_fixed.su2_ML_DE_tol0.2.txt","OneraM6_fixed.su2_ML_DE_tol0.1.txt"]
files = ["OneraM6_fixed.su2_ML_DE_tol0.1_r1.txt","OneraM6_fixed.su2_ML_DE_tol0.1_r2.txt","OneraM6_fixed.su2_ML_DE_tol0.1_r5.txt","OneraM6_fixed.su2_ML_DE_tol0.1_r10.txt"]


tol = 1e-3
#%%
#plt.close("all")x
convPlot(files, tol)
#%%
#plotCompTime(files)
#plt.show()
