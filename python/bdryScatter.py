
import matplotlib.pyplot as plt
from funs import getPlotData
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
class bdryScatterFuns:
    def scatter(v,bdryPnts,style):
        plt.scatter(v[bdryPnts][:,0],v[bdryPnts][:,1], marker= style,s = 70)
        return

    def bdryScatter(fname):
        [_,v,_,bdryPnts,_,_,_] = getPlotData(fname)
        bdryPnts = np.unique(bdryPnts)
        plt.figure()
        ax = plt.gca()
        bdryScatterFuns.scatter(v,bdryPnts,"1")    
        ax.set_aspect('equal')
        return

    def bdryScatter2(v1,bdryPnts1, fname):
        [_,v2,_,bdryPnts2,_,_,_] = getPlotData(fname)
#        bdryPnts3 = [322, 327, 337, 312]
#        bdryPnts3 = [351, 6, 19, 0, 25, 669, 656, 675, 650, 323, 352, 327, 348, 364, 389, 519, 494, 325, 673, 349, 322, 326, 350, 353, 207, 663, 156, 324, 623, 13, 260]
#        bdryPnts3 = [351, 6, 19, 25, 0, 656, 669, 675, 650, 327, 348, 285, 323, 442, 156, 493, 12]
        bdryPnts3 = []
        bdryPnts1 = np.unique(bdryPnts1)
        bdryPnts2 = np.unique(bdryPnts2)
        labels = ['Initial', 'Deformed']
        plt.figure()
        ax = plt.gca()    
        bdryScatterFuns.scatter(v1,bdryPnts1, ".")
        bdryScatterFuns.scatter(v2,bdryPnts2,".")
#        bdryScatterFuns.scatter(v2,bdryPnts3,"3")

        ax.legend(labels, bbox_to_anchor=(1.04,1), loc='upper left')
        ax.set_aspect('equal')
#        ax.set_xlim([-0.05,1.05])
#        ax.set_ylim([-0.05,1.05])
        
        return
    
    def bdryScatter3D(v, bdryPnts, markerTags, nElemsMarker, plotTag, FFD_pnts):
        fig = plt.figure()
        ax = fig.add_subplot(projection='3d')
        for i in range(len(plotTag)):
            if(plotTag[i] == "FFD"):
                ax.scatter(FFD_pnts[:,0],FFD_pnts[:,1],FFD_pnts[:,2],marker = "1")
            else:
                idx = markerTags.index(plotTag[i])        
                
                bdryPntsPlot =  np.unique(bdryPnts[sum(nElemsMarker[0:idx]):sum(nElemsMarker[0:idx+1]),:])
            
                ax.scatter(v[bdryPntsPlot][:,0],v[bdryPntsPlot][:,1],v[bdryPntsPlot][:,2],marker = "1")
#        ax.scatter(v[pnts][:,0],v[pnts][:,1],v[pnts][:,2],marker = "2", color = 'red')
#        bdryPnts = [490,499,400,409,4,904,903,949,950,459,493,404,450,469,14,948,330]
#        ax.scatter(v[bdryPnts][:,0],v[bdryPnts][:,1],v[bdryPnts][:,2],marker = "2", color = "red")
        ax.legend(plotTag)
        scaling = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
        ax.auto_scale_xyz(*[[np.min(scaling), np.max(scaling)]]*3)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        
        
