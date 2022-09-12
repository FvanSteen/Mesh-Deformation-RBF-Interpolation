# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 11:59:02 2022

@author: floyd
"""

import os
import numpy as np
import matplotlib.pyplot as plt
os.chdir('../Classical RBF' )
x = os.getcwd()

print(x)

fileName = '\square.su2'


        
fileObj = open(x+fileName, "r") 
lines = fileObj.read().splitlines()[20:29] 
fileObj.close()

arr = np.zeros((9,2))
for i in range(len(lines)):
    line = lines[i].split("\t")
    for ii in range(2):
        arr[i,ii] = float(line[ii])
        
        

plt.figure()
plt.scatter(arr[:,0],arr[:,1])
        
        
        
        
        

        


#a = lines[20]
#a = a.strip()
#a = a.split("\t")
#
#
#
#print(a)
        