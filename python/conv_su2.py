# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 00:25:33 2023

@author: floyd
"""
import matplotlib.pyplot as plt
from numpy import genfromtxt
import os
my_data = genfromtxt('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\convHist\\su2_conv.csv', delimiter=';',encoding="utf8")


#import csv
#f = open('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\MeshDeformationTool\\convHist\\su2_conv.csv', 'rt')
#reader = csv.reader(f)
#data = [] 
#for r in reader: 
#    data.append(r)
#f.close()
figPath = os.path.dirname(os.path.abspath(__file__)) + "/figs/"
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 16
plt.rcParams["figure.figsize"] = [6.4, 4.8]
plt.figure()
plt.plot(my_data[:,0],my_data[:,1], linewidth = 2, alpha = 0.6)
plt.plot(my_data[:,3],my_data[:,4], linewidth = 2, alpha = 0.6)
plt.plot(my_data[:,6],my_data[:,7], linewidth = 2, alpha = 0.6)
plt.legend(["Elastic","Pseudo","Direct"], loc = "lower right")
plt.ylabel("Lift coefficient $C_L$ [-]")
plt.xlabel("Time [min]")
plt.xlim([0,40])
plt.savefig(figPath + "Cl_conv_hd.png", dpi=400,bbox_inches='tight')
plt.figure()
plt.plot(my_data[:,0],my_data[:,2], linewidth = 2, alpha = 0.6)
plt.plot(my_data[:,3],my_data[:,5], linewidth = 2, alpha = 0.6)
plt.plot(my_data[:,6],my_data[:,8], linewidth = 2, alpha = 0.6)

plt.legend(["Elastic","Pseudo","Direct"], loc = "lower right")
plt.ylabel("Drag coefficient $C_D$ [-]")
plt.xlabel("Time [min]")
plt.xlim([0,40])
plt.savefig(figPath + "Cd_conv_hd.png", dpi=400,bbox_inches='tight')
