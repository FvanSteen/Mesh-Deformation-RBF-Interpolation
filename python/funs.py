import os
import numpy as np

def getMeshQual(fname):
    filePath = os.path.dirname(os.path.abspath(__file__))
    qFile = filePath[0:filePath.rfind('\\')+1] + "MeshDeformationTool\\meshQuals" + fname[0:-4] + "_qual.txt"
    meshQ = np.genfromtxt(qFile)
    
    print("\nMesh quality parameters of", str(fname[1:-4]) + ":" )
    print("Min mesh quality: \t", round(np.min(meshQ),5)) 
    print("Max mesh quality: \t", round(np.max(meshQ),5)) 
    print("Mean mesh quality: \t", round(np.mean(meshQ),5), "\n")
    
    
    return meshQ

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
    f_ss = np.sqrt(f_size)*f_skew
    
    return f_ss

def getMeshQuals3D(faces,vertices, alphas_0, elemType):
    
    [alphas,lambda_11, lambda_22, lambda_33] = getMeshQualParams3D(faces,vertices,elemType)
    
    tau = np.sum(alphas,axis=1)/np.sum(alphas_0, axis=1)
    f_size = np.amin(np.array([tau,1/tau]),axis=0)
    

    f_skew = np.empty(np.size(elemType))
        
    f_skew = 8/np.sum(np.sqrt(lambda_11*lambda_22*lambda_33)/alphas, axis=1)
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
        elif elemType[i] == 12:
            size = 8
            
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


def getMeshQualParams3D(faces,vertices,elemType):
    
    alphas = np.zeros(faces.shape)
    lambda_11 = np.zeros(faces.shape)
    lambda_22 = np.zeros(faces.shape)
    lambda_33 = np.zeros(faces.shape)
    
    kp1 = [1,2,3,0,7,4,5,6]
    kp3 = [3,0,1,2,5,6,7,4]
    kp4 = [4,5,6,7,0,1,2,3]
    
    # for each elemeent
    for i in range(0,faces.shape[0]):    
        
        # amount of vertices
        if elemType[i] == 12:
            size = 8
        
        # the coordinates of the element
        elem = vertices[faces[i,0:size+1]]

        for k in range(0,size):
            A = np.array([[elem[kp1[k],0]-elem[k,0], elem[kp3[k],0]-elem[k,0], elem[kp4[k],0]-elem[k,0]],
                          [elem[kp1[k],1]-elem[k,1], elem[kp3[k],1]-elem[k,1], elem[kp4[k],1]-elem[k,1]],
                          [elem[kp1[k],2]-elem[k,2], elem[kp3[k],2]-elem[k,2], elem[kp4[k],2]-elem[k,2]]])
    
            alphas[i,k] = np.linalg.det(A)
    
            tensor = np.matmul(np.transpose(A),A)
            
            lambda_11[i,k] = tensor[0,0]
            lambda_22[i,k] = tensor[1,1]
            lambda_33[i,k] = tensor[2,2]
        

    return [alphas,lambda_11,lambda_22,lambda_33]




def getPlotData(fileName):
    x = os.getcwd()
    fileObj = open(x+fileName, "r") 
    lines = fileObj.read().splitlines()

    fileObj.close()

    nMarkIndices = []
    nElemsMarks = []
    markerTags = []
    idx = 0
    FFDCtrlPntsIdx = -1
    for line in lines:
        
        if line.strip().startswith('NDIME='):
            nDim = int(line[7:].split()[0])
            nDimIdx = idx
        elif line.strip().startswith('NELEM='):
            nElem = int(line[7:].split()[0])
            nElemIdx = idx
        elif line.strip().startswith('NPOIN='):
            nPnt = int(line[7:].split()[0])
            nPntIdx = idx            
        elif line.strip().startswith('MARKER_TAG= '):
            if(line.strip().find("SEND_RECEIVE") != -1):
                save = False
            else:
                save = True
                markerTags.append(str(line[11:].split()[0]))
            
        elif line.strip().startswith('MARKER_ELEMS='):
            if(save):
                nMarkIndices.append(idx)
                nElemsMarks.append(int(line[14:].split()[0]))
        elif line.strip().startswith('FFD_CONTROL_POINTS='):
            nFFDCtrlPnts = int(line[19:].split()[0])
            FFDCtrlPntsIdx = idx+1
            
            
#        elif line.strip().startswith('MARKER_TAG= '+intBdryTag):
#            intBdryIdx = idx+1
        
        idx += 1
    print("Dimensions:\t", nDim, "\nElements:\t",nElem, "\nPoints:\t\t", nPnt)

    
    nNodesElem = 8
    
    f = np.empty((nElem,nNodesElem), dtype = int)
    f[:] = -1

    elemType = np.empty(nElem)
    
    maxNodesElem = 3

    for i in range(nElem):
        lineData = lines[i+nElemIdx+1].strip().split()
        
        if lineData[0] == '5':
            f[i,0:3] = lineData[1:4]
            elemType[i] = 5
            
        elif lineData[0] == '9':
            f[i,0:4] = lineData[1:5]
            elemType[i] = 9
            if maxNodesElem < 4:
                maxNodesElem = 4
        elif lineData[0] == "12":
            f[i,0:8] = lineData[1:9]
            elemType[i] = 12
            if maxNodesElem < 8:
                maxNodesElem = 8
        elif lineData[0] == "10":
            f[i,0:4] = lineData[1:5]
            elemType[i] = 10
            if maxNodesElem < 4:
                maxNodesElem = 4
#        f[i,:] = lineData[1:nNodesElem+1]
#        print(f[i,:])
    
   
    f = np.delete(f, np.s_[maxNodesElem:nNodesElem], axis=1)  
   
    
        
    
    v = np.empty((nPnt,nDim))
    for i in range(nPnt):        
        lineData = lines[i+nPntIdx+1].strip().split()
        v[i,:] = lineData[:nDim]

            

    bdryType = int(lines[nMarkIndices[0]+1].strip().split()[0])
    
    if bdryType == 5:
        cols = 3
    elif bdryType == 9:
        cols = 4
    elif bdryType == 3:
        cols = 2
        
    bdry = np.empty((np.sum(nElemsMarks),cols), dtype = int)
    for ii in range(len(nMarkIndices)):        
        for j in range(nElemsMarks[ii]):
#            # first strip and then append
            
            lineData = lines[nMarkIndices[ii]+j+1].strip().split()
            
            bdry[int(np.sum(nElemsMarks[0:ii])+j),:] = np.array(lineData[1:], dtype = int)        
    
#    bdry = np.unique(bdry)
    
    
    if FFDCtrlPntsIdx != -1:
        FFD_pnts = np.empty((nFFDCtrlPnts,nDim),dtype = float)
        for i in range(FFDCtrlPntsIdx, FFDCtrlPntsIdx+nFFDCtrlPnts):
            lineData = lines[i].strip().split()
            FFD_pnts[i-FFDCtrlPntsIdx] = np.array(lineData[3:],dtype = float)
            
    else:
        FFD_pnts = np.empty((0,0))
    
            
    return [f,v,elemType,bdry, markerTags, nElemsMarks, FFD_pnts]