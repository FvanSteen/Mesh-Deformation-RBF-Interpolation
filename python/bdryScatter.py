# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 20:41:37 2022

@author: floyd
"""

import matplotlib.pyplot as plt
from funs import getPlotData
from mpl_toolkits.mplot3d import Axes3D
class bdryScatterFuns:
    def scatter(v,bdryPnts,style):
        plt.scatter(v[bdryPnts][:,0],v[bdryPnts][:,1], marker= style)
        return

    def bdryScatter(fname):
        [_,v,_,bdryPnts] = getPlotData(fname)
        plt.figure()
        ax = plt.gca()
        bdryScatterFuns.scatter(v,bdryPnts,"1")    
        ax.set_aspect('equal')
        return

    def bdryScatter2(v1,bdryPnts1, fname):
        [_,v2,_,bdryPnts2] = getPlotData(fname)
#        bdryPnts3 = [14245,15478,15079,15469,14458,15078,14891,15439,15396,14424,14350,14873,15476,15284]
#        bdryPnts3 = [15468]
#        bdryPnts3 = [322,10,348,650,654,675,0,353,390,337,208]
#        bdryPnts3 = [97008,97662,98319, 2587,    2, 2400,    6, 2129,    3,97385, 2487, 2013,97818, 1550,97994,    5, 1383, 2347, 1750,97519,97907, 2707, 2329,98095,97758,97431, 2293,97864, 1402, 2374, 1848,97475,97230, 1662,98205,97304,97618, 1617,97111, 1793, 1346,97841, 2361,97940,98041, 1440, 1511, 2524,97150, 2070, 2312,97567,97346, 1707, 2201, 1410, 1579,97964,97718, 2444,97059, 1769,97881, 2337,97267,97408,97700,97646,98274, 1639, 2269, 2673,97455,97497,97788, 1932,97185,97894,97083,97923, 1818, 1331, 1729, 1364, 2387,97680,98149,97871, 2231,97855,98240, 2549, 2466, 2505, 1871, 1414, 1471, 2623, 2693, 1497, 1527,97657, 1565, 2368,97325,97596, 1781,97738, 1684,97830, 2355, 2424,97132,97095,98296,98017,97035]
#        bdryPnts3 = [14424]
        bdryPnts3 = [326,348,155,520,672,0,441,353]
        
        labels = ['Initial Mesh', 'Deformed Mesh']
        plt.figure()
        ax = plt.gca()    
        bdryScatterFuns.scatter(v1,bdryPnts1,"1")
        bdryScatterFuns.scatter(v2,bdryPnts2,"2")
        bdryScatterFuns.scatter(v2,bdryPnts3,"3")
#        plt.scatter( 0.0488328, -0.0760844)
#        plt.scatter(0.0548919, -0.0884848)
#        plt.scatter(0.0482768, -0.0759677)
#        plt.scatter(0.0571281, -0.0778245)
#        plt.scatter(0.0499529, -0.0738999)
#        plt.scatter(0.0498472, -0.0741069)
        plt.scatter(0.0488328, -0.0760844)
        plt.scatter(0.0548919, -0.0884848)
        plt.scatter(0.0571281, -0.0778245)
        plt.scatter(0.0498472, -0.0741069)
        ax.legend(labels)
        ax.set_aspect('equal')
        return
    
    def bdryScatter3D(fname):
        [_,v,_,bdryPnts] = getPlotData(fname)
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        ax.scatter(v[bdryPnts][:,0],v[bdryPnts][:,1],v[bdryPnts][:,2],marker = "1")
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        
        
