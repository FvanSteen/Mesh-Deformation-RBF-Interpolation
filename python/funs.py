import os
import numpy as np


def getMeshQuals(faces,vertices, alphas_0, elemType):
    
    [alphas,lambda_11, lambda_22,lambda_12] = getMeshQualParams(faces,vertices,elemType)
    
    tau = np.sum(alphas,axis=1)/np.sum(alphas_0, axis=1)
    f_size = np.amin(np.array([tau,1/tau]),axis=0)
    quadIdx = np.where(elemType == 9)
    if np.size(quadIdx) != 0:
        startQuadIdx = np.where(elemType == 9)[0][0]
    else:
        startQuadIdx = np.size(elemType)
    
    f_skew = np.empty(np.size(elemType))
    f_skew[0:startQuadIdx] = np.sqrt(3)*alphas[0:startQuadIdx,0]/(lambda_11[0:startQuadIdx,0]+lambda_22[0:startQuadIdx,0]-lambda_12[0:startQuadIdx,0])
    f_skew[startQuadIdx:] = np.abs(4/np.sum(np.sqrt(lambda_11[startQuadIdx:]*lambda_22[startQuadIdx:])/alphas[startQuadIdx:], axis=1))
#    if elementType == 5:
#        f_skew = np.sqrt(3)*alphas[:,0]/(lambda_11[:,0]+lambda_22[:,0]-lambda_12[:,0])
#        
#    elif elementType == 9:
#        f_skew = 4/np.sum(np.sqrt(lambda_11*lambda_22)/alphas, axis=1)
    f_ss = np.sqrt(f_size)*f_skew
    return f_ss
    
def getMeshQualParams(faces, vertices, elemType):
    alphas = np.zeros(faces.shape)
    lambda_11 = np.zeros(faces.shape)
    lambda_22 = np.zeros(faces.shape)
    lambda_12 = np.zeros(faces.shape)
    
    for i in range(0,faces.shape[0]):    
        if elemType[i] == 5:
            size = 3
        elif elemType[i] == 9:
            size = 4
            
        elem = vertices[faces[i,0:size+1]]
        
        for ii in range(0,size):
            iip1 = (ii+1)%(size)
            iip3= (ii+(size-1))%(size)
            
            A = np.array([[elem[iip1,0]-elem[ii,0], elem[iip3,0]-elem[ii,0]], [elem[iip1,1]-elem[ii,1], elem[iip3,1]-elem[ii,1]]])

            alphas[i,ii] = np.linalg.det(A)
            
            tensor = np.matmul(np.transpose(A),A)
            
            lambda_11[i,ii] = tensor[0,0]
            lambda_22[i,ii] = tensor[1,1]
            lambda_12[i,ii] = tensor[0,1]
                      
            
    return [alphas,lambda_11,lambda_22,lambda_12]

def getPlotData(fileName):
    x = os.getcwd()
    fileObj = open(x+fileName, "r") 

    lines = fileObj.read().splitlines()
    fileObj.close()

    nMarkIndices = []
    nMarks = []
    idx = 0
    for line in lines:
        
        if line.strip().startswith('NDIME= '):
            nDim = int(line[7:])
            nDimIdx = idx
        elif line.strip().startswith('NELEM= '):
            nElem = int(line[7:])
            nElemIdx = idx
        elif line.strip().startswith('NPOIN= '):
            nPnt = int(line[7:])
            nPntIdx = idx
        elif line.strip().startswith('MARKER_ELEMS= '):
            nMarkIndices.append(idx)
            nMarks.append(int(line[14:]))
            
#        elif line.strip().startswith('MARKER_TAG= '+intBdryTag):
#            intBdryIdx = idx+1
        
        idx += 1
    
    
    print("Dimensions:\t", nDim, "\nElements:\t",nElem, "\nPoints:\t\t", nPnt)
#    print(nDimIdx,nElemIdx,nPntIdx)
    
#    nLine = lines[intBdryIdx]
#    nInElems = int(nLine[14:])
#    f_in = np.empty([1,2*nInElems],dtype=int)
#    
#    for i in range(nInElems):
#        lineData = lines[intBdryIdx+1+i].strip().split('\t')    
#        f_in[0,2*i:2*i+2] = lineData[1:3]
    
    # Change the 3 and 4 here to be adjustable to the type of elements used in the mesh
    
#    if elementType == 5:
#        nNodesElem = 3
#    elif elementType == 9:
#        nNodesElem = 4
#    elif elementType == 12:
#        nNodesElem = 8
    
    nNodesElem = 4
    
    f = np.empty((nElem,nNodesElem),dtype=int)
    f[:] = -1
    elemType = np.empty(nElem)
    
    
#    print(f[0,:])
    for i in range(nElem):
        lineData = lines[i+nElemIdx+1].strip().split()
#        print(lineData)
#        print(lineData[1:nNodesElem+1])
#        print(lineData[1:4])
#        print(f[0,:])
#        print(lineData[0])
        if lineData[0] == '5':
            f[i,0:3] = lineData[1:4]
            elemType[i] = 5
        elif lineData[0] == '9':
            f[i,0:4] = lineData[1:5]
            elemType[i] = 9
#        f[i,:] = lineData[1:nNodesElem+1]
#        print(f[i,:])
        
    
        
    
    v = np.empty((nPnt,nDim))
    for i in range(nPnt):
        lineData = lines[i+nPntIdx+1].strip().split()
        v[i,:] = lineData[:nDim]

    bdry = np.array([],dtype=int)        
    for ii in range(len(nMarkIndices)):
        for j in range(nMarks[ii]):
            # first strip and then append
            lineData = lines[nMarkIndices[ii]+j+1].strip().split()
            
            bdry = np.append(bdry,np.array(lineData[1:],dtype=int))
#            print(lines[nMarkIndices[ii]+j+1])
    bdry = np.unique(bdry)
    
            
    return [f,v,elemType,bdry]