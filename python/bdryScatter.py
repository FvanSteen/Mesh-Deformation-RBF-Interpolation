# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 20:41:37 2022

@author: floyd
"""

import matplotlib.pyplot as plt
from funs import getPlotData

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
        bdryPnts3 = [15478]
        plt.figure()
        ax = plt.gca()    
        bdryScatterFuns.scatter(v1,bdryPnts1,"1")
        bdryScatterFuns.scatter(v2,bdryPnts2,"2")
        bdryScatterFuns.scatter(v2,bdryPnts3,"3")
        ax.set_aspect('equal')
        return
