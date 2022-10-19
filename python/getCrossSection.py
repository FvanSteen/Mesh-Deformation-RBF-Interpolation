import numpy as np
import math

def getCrossSection(f,v,z_cut):
    ver = []
    
    # for loop for each element
    #450,675
    for elems in range(len(f)):    
        if(elems%100 == 0):
            print(elems)
        # verts is an array containing subarrays with the coordinates of each face of the element
        verts = [np.array([v[f][elems][0],v[f][elems][1],v[f][elems][5],v[f][elems][4]]),
              np.array([v[f][elems][3],v[f][elems][2],v[f][elems][6],v[f][elems][7]]),
              np.array([v[f][elems][3],v[f][elems][0],v[f][elems][4],v[f][elems][7]]),
              np.array([v[f][elems][2],v[f][elems][1],v[f][elems][5],v[f][elems][6]]),
              np.array([v[f][elems][0],v[f][elems][1],v[f][elems][2],v[f][elems][3]]),
              np.array([v[f][elems][4],v[f][elems][5],v[f][elems][6],v[f][elems][7]])]
        
        
        # Lists for storing x and y-coordinate
        xcoords = [];
        ycoords = [];
        
        # for each fase determine wheter it crosses with z_cut
        for i in range(len(verts)):    
            if np.all(verts[i][:,2] == z_cut):
                ver.append(verts[i][:,0:2])
                        
            if np.any(verts[i][:,2] < z_cut) and np.any(verts[i][:,2] > z_cut): # if on both sides of cutline
                
                for ii in range(len(verts[i])):      # for edges
                    if ii < len(verts[i])-1 and  np.any(verts[i][ii:ii+2,2] < z_cut) and np.any(verts[i][ii:ii+2,2] > z_cut):
#                        print(ii, verts[i][ii,:], verts[i][ii+1,:])
                        if verts[i][ii+1,2]-verts[i][ii,2] < 0:
                            fact = verts[i][ii,2]-z_cut
                        else:
                            fact = z_cut-verts[i][ii,2]
#                        print(fact/(verts[i][ii+1,2]-verts[i][ii,2]))
                        x = (verts[i][ii+1,0]-verts[i][ii,0])*fact/abs(verts[i][ii+1,2]-verts[i][ii,2]) + verts[i][ii,0]
                        y = (verts[i][ii+1,1]-verts[i][ii,1])*fact/abs(verts[i][ii+1,2]-verts[i][ii,2]) + verts[i][ii,1]
#                        print(fact, fact/abs(verts[i][ii+1,2]-verts[i][ii,2]),x,y)
                        
                        xcoords.append(x)
                        ycoords.append(y)
                    elif ii== len(verts[i])-1 and np.any(np.array([verts[i][ii,2], verts[i][0,2]]) <z_cut) and np.any(np.array([verts[i][ii,2], verts[i][0,2]]) > z_cut):
#                        print(ii, verts[i][0,:], verts[i][ii,:])
                        if verts[i][ii,2]-verts[i][0,2] < 0:
#                            fact = z_cut-verts[i][0,2]
                            fact = verts[i][0,2] - z_cut
                        else:
#                            fact = verts[i][ii,2]-z_cut
                            fact = z_cut-verts[i][0,2]
                            
#                        print(fact, fact/abs(verts[i][ii,2]-verts[i][0,2]),x,y)
                        x = (verts[i][ii,0]-verts[i][0,0])*fact/abs(verts[i][ii,2]-verts[i][0,2]) + verts[i][0,0]
                        y = (verts[i][ii,1]-verts[i][0,1])*fact/abs(verts[i][ii,2]-verts[i][0,2]) + verts[i][0,1]
                        
                        
                        xcoords.append(x)
                        ycoords.append(y)
    
        
        
        if len(xcoords) > 1:
            
            coords = np.array([xcoords,ycoords])
             
            coords = np.unique(np.transpose(coords),axis=0)
            
            midpoint = np.sum(coords,axis=0)/np.size(coords,axis=0)
            vecs = coords-midpoint
            theta = np.zeros((len(vecs),1))
            for i in range(len(coords)):
                theta[i] = math.atan2(vecs[i,1],vecs[i,0])*180/np.pi
            
            coords = np.column_stack((coords,theta))
            idx = np.argsort(coords[:,2])
            coords = coords[idx]
            
            ver.append(coords[:,0:2])
    return ver
            