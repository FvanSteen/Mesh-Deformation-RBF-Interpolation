import numpy as np
import math

def getCrossSection(f,v,z_cut,dim,meshQual,f_in):
    ver = []
    mq= []
    dims = [0,1,2]
    dims.remove(dim)
    ver_in = np.array([])
    # for loop for each element
    #450,675
    for elems in range(len(f)):    
        
        if set(f[elems]).issubset(set(f_in)) and np.any(v[f][elems] < z_cut) and np.any(v[f][elems] > z_cut):
            ver_in = np.array([v[f][elems][0:4,0:2]])
        
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
        coord1 = [];
        coord2 = [];
        
        # for each fase determine wheter it crosses with z_cut
        for i in range(len(verts)):    
            if np.all(verts[i][:,dim] == z_cut):
                ver.append(verts[i][:,0:2])
                mq.append(meshQual[elems])
                        
            if np.any(verts[i][:,dim] < z_cut) and np.any(verts[i][:,dim] > z_cut): # if on both sides of cutline
                
                for ii in range(len(verts[i])):      # for edges
                    if ii < len(verts[i])-1 and  np.any(verts[i][ii:ii+2,dim] < z_cut) and np.any(verts[i][ii:ii+2,dim] > z_cut):
#                        print(ii, verts[i][ii,:], verts[i][ii+1,:])
                        if verts[i][ii+1,dim]-verts[i][ii,dim] < 0:
                            fact = verts[i][ii,dim]-z_cut
                        else:
                            fact = z_cut-verts[i][ii,dim]
#                        print(fact/(verts[i][ii+1,2]-verts[i][ii,2]))
                        c1 = (verts[i][ii+1,dims[0]]-verts[i][ii,dims[0]])*fact/abs(verts[i][ii+1,dim]-verts[i][ii,dim]) + verts[i][ii,dims[0]]
                        c2 = (verts[i][ii+1,dims[1]]-verts[i][ii,dims[1]])*fact/abs(verts[i][ii+1,dim]-verts[i][ii,dim]) + verts[i][ii,dims[1]]
#                        print(fact, fact/abs(verts[i][ii+1,2]-verts[i][ii,2]),x,y)
                        
                        coord1.append(c1)
                        coord2.append(c2)
                    elif ii== len(verts[i])-1 and np.any(np.array([verts[i][ii,dim], verts[i][0,dim]]) <z_cut) and np.any(np.array([verts[i][ii,dim], verts[i][0,dim]]) > z_cut):
#                        print(ii, verts[i][0,:], verts[i][ii,:])
                        if verts[i][ii,dim]-verts[i][0,dim] < 0:
#                            fact = z_cut-verts[i][0,2]
                            fact = verts[i][0,dim] - z_cut
                        else:
#                            fact = verts[i][ii,2]-z_cut
                            fact = z_cut-verts[i][0,dim]
                            
#                        print(fact, fact/abs(verts[i][ii,2]-verts[i][0,2]),x,y)
                        c1 = (verts[i][ii,dims[0]]-verts[i][0,dims[0]])*fact/abs(verts[i][ii,dim]-verts[i][0,dim]) + verts[i][0,dims[0]]
                        c2 = (verts[i][ii,dims[1]]-verts[i][0,dims[1]])*fact/abs(verts[i][ii,dim]-verts[i][0,dim]) + verts[i][0,dims[1]]
                        
                        
                        coord1.append(c1)
                        coord2.append(c2)
    
        
        ## todo change the way coords are assigned here. 
        if len(coord1) > 1:
            
            coords = np.array([coord1,coord2])
             
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
            mq.append(meshQual[elems])
          
    return [ver,ver_in,mq]
            