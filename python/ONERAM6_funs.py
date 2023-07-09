import os
import numpy as np
import matplotlib.pyplot as plt

def rewriteSU2File(fName):
    
    x = os.getcwd()
        
    fileObj = open(x+fName, "r") 
    output = open(x+"\\OneraM6_fixed.su2", "w")
    
    lines = fileObj.read().splitlines()
    pnts = False
    cnt = 0
    for line in lines:
        if line.strip().startswith('NPOIN= '):
            nPnts = int(line[7:].split()[0])
            pnts = True
            output.write('NPOIN= '+ str(nPnts)  + "\n")
        elif pnts and cnt < nPnts:
            output.write("\t".join(line.split()[0:-1]) + "\n" )
            cnt+= 1
        else:
            output.write(line + "\n")
    output.close()
    fileObj.close()
    print("done rewriting file in right format") 

def getFFDdef(FFD_pnts):
    ymin = np.min(FFD_pnts[:,1])
    ymax = np.max(FFD_pnts[:,1])    
    
    root_pnts = FFD_pnts[np.where(FFD_pnts[:,1]==ymin),:][0]
    tip_pnts = FFD_pnts[np.where(FFD_pnts[:,1]==ymax),:][0]
    
    root_chord = np.max(root_pnts[:,0]) - np.min(root_pnts[:,0])
    tip_chord = np.max(tip_pnts[:,0]) - np.min(tip_pnts[:,0])    
    
    taper = tip_chord/root_chord
    mac = 2/3*root_chord*((1+taper+taper**2)/(1+taper))
    print("MAC",mac)
    TE_loc_root = np.max(root_pnts[:,0])
    TE_loc_tip = np.max(tip_pnts[:,0])
        
#    dz_wing = mac*np.sin(FFD_pnts[:,1]/(2*ymax)*np.pi)
    
    disp = np.zeros(np.shape(FFD_pnts), dtype = float)
    
    dz_wing = 0.25*(2*ymax)*FFD_pnts[:,1]*(1-np.cos(FFD_pnts[:,1]*np.pi)/(2*ymax))
    
#    phi_0 = 30/180*np.pi
#    
#    for i in range(np.shape(FFD_pnts)[0]):    
#        x_te = TE_loc_root + FFD_pnts[i,1]/(2*ymax) *(TE_loc_tip-TE_loc_root)
#        dx = FFD_pnts[i,0] - x_te
#        dz = FFD_pnts[i,2]
#        
#        phi = phi_0*np.sin(FFD_pnts[i,1]/(2*ymax) * np.pi)
#        
#        xnew = dx*np.cos(phi)-dz*np.sin(phi) + x_te
#        znew = dx*np.sin(phi)+dz*np.cos(phi)
#        disp[i,0] = xnew - FFD_pnts[i,0]
#        disp[i,2] = znew - FFD_pnts[i,2]

    disp[:,2] += dz_wing
    
    return disp


def getOutputSU2File(FFD_pnts, v_init, bdryPnts_init, markerTags, nElemsMarks, wingTags):
    
    x = os.getcwd()
    o = open(x + "\\OneraM6_FFD.su2", "w")
    
    xElem = 10
    yElem = 8
    zElem = 1
    
    xNodes = xElem +1
    yNodes = yElem +1
    zNodes = zElem +1
    idx = []
    for i in wingTags:
        idx.append(markerTags.index(i))
    
    
    wingPnts = np.empty((0,3), dtype = int)
    for i in range(len(idx)):
        wingPnts = np.vstack((wingPnts, bdryPnts_init[sum(nElemsMarks[0:idx[i]]):sum(nElemsMarks[0:idx[i]+1])]))      
    
    
    wingNodes = np.unique(wingPnts)
    
    o.write("NDIME= 3 \n")
    o.write("NELEM= 0\n")# + str(xElem*yElem*zElem) + "\n")
    o.write("NPOIN= " +  str(xNodes*yNodes*zNodes + np.size(wingNodes)) + "\n")
#    for i in range(np.shape(FFD_pnts)[0]):
#        o.write("\t".join(FFD_pnts) + "\n" )
    
    
    
    for c1,c2,c3 in FFD_pnts:
        o.write("{0}\t{1}\t{2}".format(c1,c2,c3) + "\n")
        
    for c1,c2,c3 in v_init[wingNodes]:
        o.write("{0}\t{1}\t{2}".format(c1,c2,c3) + "\n")
    
    
    o.write("NMARK= 6\n")
    o.write("MARKER_TAG= FRONT\n")
    o.write("MARKER_ELEMS= "+str(zElem*yElem) + "\n")
    
    for i in range(yElem):
        o.write( "9\t" +  str(i*zNodes) + "\t" +  str((i*zNodes)+1)  + "\t" +  str((i+1)*zNodes+1) + "\t" +  str((i+1)*zNodes) + "\n")
        
    o.write("MARKER_TAG= BACK\n")
    o.write("MARKER_ELEMS= "+str(zElem*yElem) + "\n")
    
    for i in range(yElem):
        o.write("9\t"+ str(i*zNodes+ zNodes*yNodes*(xNodes-1)) + "\t" +  str((i*zNodes)+1 + zNodes*yNodes*(xNodes-1)) + "\t" +  str((i+1)*zNodes+1 +zNodes*yNodes*(xNodes-1)) + "\t" +  str((i+1)*zNodes + zNodes*yNodes*(xNodes-1)) +  "\n")
    
    
    
    o.write("MARKER_TAG= ROOT\n")
    o.write("MARKER_ELEMS= "+str(zElem*xElem) + "\n")
    
    for i in range(xElem):
        o.write("9\t"+ str(i*yNodes*zNodes)+ "\t" + str((i+1)*yNodes*zNodes)+ "\t" +  str((i+1)*yNodes*zNodes+1) + "\t" + str((i)*yNodes*zNodes+1) + "\n" )
   
    o.write("MARKER_TAG= TIP\n")
    o.write("MARKER_ELEMS= "+str(zElem*xElem) + "\n")
    
    for i in range(xElem):
        o.write( "9\t" + str(i*yNodes*zNodes + zNodes*(yNodes-1)) + "\t" + str((i+1)*yNodes*zNodes + zNodes*(yNodes-1)) + "\t" + str((i+1)*yNodes*zNodes+1 + zNodes*(yNodes-1)) + "\t" + str((i)*yNodes*zNodes+1 + zNodes*(yNodes-1))  + "\n")
    
    # lower
    o.write("MARKER_TAG= LOWER\n")
    o.write("MARKER_ELEMS= "+str(yElem*xElem) + "\n")
    for i in range(xElem):
        for ii in range(yElem):
            o.write("9\t"+ str(ii*zNodes + i*yNodes*zNodes) + "\t" +  str(ii*zNodes + (i+1)*yNodes*zNodes) + "\t" +  str((ii+1)*zNodes + (i+1)*yNodes*zNodes) + "\t" + str((ii+1)*zNodes + i*yNodes*zNodes) + "\n")
    #   
    o.write("MARKER_TAG= UPPER\n")
    o.write("MARKER_ELEMS= "+str(yElem*xElem) + "\n")
    # upper     
    for i in range(xElem):
        for ii in range(yElem):
            o.write("9\t" + str(ii*zNodes + i*yNodes*zNodes + 1)  + "\t" +  str(ii*zNodes + (i+1)*yNodes*zNodes+ 1)  + "\t" +  str((ii+1)*zNodes + (i+1)*yNodes*zNodes+ 1)  + "\t" + str((ii+1)*zNodes + i*yNodes*zNodes+ 1) + "\n")
    o.close()
    
def getDefFile(defFileName, index, disp):
    x = os.getcwd()
    o = open(x + defFileName, "w")
    
    
    
    for idx,c1,c2,c3 in zip(np.transpose(index),disp[:,0],disp[:,1],disp[:,2]):
        o.write("{0}\t{1}\t{2}\t{3}".format(idx, c1,c2,c3) + "\n")
    o.close()
    
def scatteredWing(v_init, v_def, bdryPnts):
    
#    startIdxWing = np.sum(nElemsMarker[:])
    bdryNodes = np.unique(bdryPnts)

    wingPnts_init = np.delete(v_init, bdryNodes, axis=0)
    wingPnts_def = np.delete(v_def, bdryNodes, axis=0)
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.scatter(wingPnts_init[:,0],wingPnts_init[:,1],wingPnts_init[:,2])
    ax.scatter(wingPnts_def[:,0],wingPnts_def[:,1],wingPnts_def[:,2])
    ax.legend(["Undeformed", "Deformed"])
    
    scaling = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    ax.auto_scale_xyz(*[[np.min(scaling), np.max(scaling)]]*3)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    
    
def getDisp(v_FFD, v_def, bdryPnts_FFD, bdryPnts_init, markerTags, wingTags, nElemsMarks): 
    bdryNodes = np.unique(bdryPnts_FFD)    
    idx = []
    for i in wingTags:
        idx.append(markerTags.index(i))
    
    wingPnts = np.empty((0,3), dtype = int)
    for i in range(len(idx)):
        wingPnts = np.vstack((wingPnts, bdryPnts_init[sum(nElemsMarks[0:idx[i]]):sum(nElemsMarks[0:idx[i]+1])]))      
    
    wingNodes = np.unique(wingPnts)
    
    disp = v_def-v_FFD
    disp = np.delete(disp, bdryNodes, axis=0)
    return [wingNodes, disp]
    
    
    
    
    

    
    