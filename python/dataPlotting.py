import os
from matplotlib import pyplot as plt
import numpy as np

os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\')

def getData(fname):
    x = os.getcwd()
    
    fileObj = open(x+"\\data\\"+fname, "r") 
 
    lines = fileObj.read().splitlines()
    
    steps = np.empty((len(lines),1))    
    time = np.empty((len(lines),1))    
    qmin = np.empty((len(lines),1))    
    
    for i in range(len(lines)):
        data = lines[i].split()
        steps[i] = data[0]
        time[i] = data[1]
        qmin[i] = data[2]
        
    return[steps,time,qmin]
        
        
#fname = ["none_none_none.txt", "ps_none_none.txt","ds_none_none.txt"]
#fname = ["ds_none_none2.txt", "ds_periodic_none2.txt", "ds_fixed_none2.txt", "ds_moving_none2.txt"]
fname = ["none_none_none2.txt", "none_periodic_none2.txt"]

plt.close("all")
labels = ["none", "periodic", "fixed", "moving"]
plt.figure()
for i in fname:
    [steps,time,qmin] = getData(i)
    plt.plot(steps,qmin, marker = '.')
    plt.xticks(range(0,21,2))
    plt.xlabel("steps [-]")
    plt.ylabel("min qual [-]")    
    plt.legend(labels)
    
    
plt.figure()
for i in fname:
    [steps,time,qmin] = getData(i)
    plt.plot(steps,time, marker = '.')
    plt.xticks(range(0,21,2))
    plt.xlabel("steps [-]")
    plt.ylabel("comp time [-]")    
    plt.legend(labels)



#labels = ["std", "ps","ds"]
#plt.figure()
#for i in fname:
#    [steps,time,qmin] = getData(i)
#    plt.plot(steps,qmin, marker = '.')
#    plt.xticks(range(0,21,2))
#    plt.xlabel("steps [-]")
#    plt.ylabel("min qual [-]")    
#    plt.legend(labels)
#    
#plt.figure()
#for i in fname:
#    [steps,time,qmin] = getData(i)
#    plt.plot(steps,time,marker = '.')
#    plt.xticks(range(0,21,2))
#    plt.xlabel("steps [-]")
#    plt.ylabel("comp time [s]")    
#    plt.legend(labels)