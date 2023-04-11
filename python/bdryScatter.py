
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
        bdryScatterFuns.scatter(v,bdryPnts,".")  
        ax.set_aspect('equal')
        return

    def bdryScatter2(v1,bdryPnts1, fname):
        [_,v2,_,bdryPnts2,_,_,_] = getPlotData(fname)
        bdryPnts1 = np.unique(bdryPnts1)
        bdryPnts2 = np.unique(bdryPnts2)
        labels = ['Initial', 'Deformed']
        plt.figure()
        ax = plt.gca()    
        bdryScatterFuns.scatter(v1,bdryPnts1, ".")
        bdryScatterFuns.scatter(v2,bdryPnts2,".")
        ax.legend(labels, bbox_to_anchor=(1.04,1), loc='upper left')
        ax.set_aspect('equal')       
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
                
        bdryPntsPlot = np.array([453,546,456,556,446,  9,543,909,444,553,990,554,300,799,390,609,199,800, 60,293,704,792,289,770, 11,260,502,659,307, 12,560,610])
        ax.scatter(v[bdryPntsPlot][:,0],v[bdryPntsPlot][:,1],v[bdryPntsPlot][:,2],marker = "1")
        ax.legend(plotTag)
#        scaling = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
#        ax.auto_scale_xyz(*[[np.min(scaling), np.max(scaling)]]*3)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.show()
        
        
