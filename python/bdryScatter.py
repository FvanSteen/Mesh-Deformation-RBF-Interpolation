
import matplotlib.pyplot as plt
from funs import getPlotData
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
class bdryScatterFuns:
    def scatter(v,bdryPnts,style):
        plt.scatter(v[bdryPnts][:,0],v[bdryPnts][:,1], marker= style,s = 70)
        return

    def bdryScatter(fname, v, bdryPnts, init):
        bdryPnts = np.unique(bdryPnts)
        plt.figure()
        ax = plt.gca()
        if(init):
            [_,v_init,_,_,_,_,_] = getPlotData(fname[:-8] + ".su2")
            bdryScatterFuns.scatter(v_init,bdryPnts,".")  
            
        bdryScatterFuns.scatter(v,bdryPnts,".")  
        
        if(init):
            ax.legend(['Initial','Deformed'], bbox_to_anchor=(1.04,1), loc='upper left')
        
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
            
                ax.scatter(v[bdryPntsPlot][:,0],v[bdryPntsPlot][:,1],v[bdryPntsPlot][:,2],marker = "1", alpha = 0.2)

        bdryPntsPlot = np.array([36480, 36482, 36484, 36486, 36488, 36490, 36492, 36494, 36496, 36498, 36500, 36502, 36504, 36506, 36508, 36510, 36512, 36514, 36516, 36518, 36520, 36522, 36524, 36526, 36528, 36530, 36532, 36534, 36536, 36538, 36540, 36542, 36544, 36546, 36548, 36550, 36552, 36554, 36556, 36558, 36560, 36562, 36564, 36566, 36568, 36570, 36572, 36574, 36576, 36578, 36580, 36582, 36584, 36586, 36588, 36590, 36592, 36594, 36596, 36598, 36600, 36602, 36604, 36606, 36608, 36610, 36612, 36614, 36616, 36618, 36620, 36622, 36624, 36626, 36628, 38684, 38685, 38686, 38687, 38688, 38689, 38690, 38691, 38692, 38693, 38694, 38695, 38696, 38697, 38698, 38699, 38700, 38701, 38702, 38703, 38704, 38705, 38706, 38707, 38708, 38709, 38710, 38711, 38712, 38713, 38714, 38715, 38716, 38717, 38718, 38719, 38720, 38721, 38722, 38723, 38724, 38725, 38726, 38727, 38728, 38729, 38730, 38731, 38732, 38733, 38734, 38735, 38736, 38737, 38738, 38739, 38740, 38741, 38742, 38743, 38744, 38745, 38746, 38747, 38748, 38749, 38750, 38751, 38752, 38753, 38754, 38755, 38756, 38757, 38758, 44880, 44884, 44888, 44892, 44896, 44900, 44904, 44908, 44912, 44916, 44920, 44924, 44928, 44932, 44936, 44940, 44944, 44948, 44952, 44956, 44960, 44964, 44968, 44972, 44976, 44980, 44984, 44988, 44992, 44996, 45000, 45004, 45008, 45012, 45016, 45020, 45024, 45028, 45032, 47200, 47202, 47204, 47206, 47208, 47210, 47212, 47214, 47216, 47218, 47220, 47222, 47224, 47226, 47228, 47230, 47232, 47234, 47236, 47238, 47240, 47242, 47244, 47246, 47248, 47250, 47252, 47254, 47256, 47258, 47260, 47262, 47264, 47266, 47268, 47270, 47272, 47274, 47276, 36630, 38759, 45036, 47278, 252510, 252566, 252570, 252599])
        ax.scatter(v[bdryPntsPlot][:,0],v[bdryPntsPlot][:,1],v[bdryPntsPlot][:,2],marker = "1", alpha =.5, color = 'green')
        bdryPntsPlot = np.array([38745, 47234, 36490, 45036, 36530, 38702, 47270,127214, 126210, 107250, 111075, 74865, 175980, 95237, 157304, 100031, 146100, 89967, 162285, 143024, 122352, 159004, 91023, 158850, 130852, 132773, 123117, 138960, 123763, 84357, 100830, 95577, 99674, 92381, 155247, 153513, 153139, 95900, 152204, 147597, 137397, 108579])
        ax.scatter(v[bdryPntsPlot][:,0],v[bdryPntsPlot][:,1],v[bdryPntsPlot][:,2], s = 50, marker = "1", alpha =1, color = 'red')
        
        bdryPntsPlot = np.array([20444])
        ax.scatter(v[bdryPntsPlot][:,0],v[bdryPntsPlot][:,1],v[bdryPntsPlot][:,2], s = 50, marker = "1", alpha =1, color = 'm')
        ax.scatter(0.3, -0.261966, 0.0831353, s = 50, marker = ".", alpha =1, color = 'blue')
        ax.scatter(0.3,-0.293191,0.0664157, s = 50, marker = ".", alpha =1, color = 'm')
        ax.scatter(0.299, -0.262649, 0.0843223, s = 50, marker = ".", alpha =1, color = 'blue')
        ax.scatter(0.3,-0.293191,0.0843223, s = 50, marker = ".", alpha =1, color = 'm')
#        bdryPntsPlot = np.array([107250])
#        ax.scatter(v[bdryPntsPlot][:,0],v[bdryPntsPlot][:,1],v[bdryPntsPlot][:,2],marker = "1", alpha = 1, color = 'red')
        
        ax.legend(plotTag)
#        scaling = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
#        ax.auto_scale_xyz(*[[np.min(scaling), np.max(scaling)]]*3)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        scaling = 1.5
        fig.set_size_inches(scaling*6.4, scaling*4.8)
        plt.show()
        
        
