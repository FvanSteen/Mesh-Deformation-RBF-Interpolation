import os
from bdryScatter import bdryScatterFuns
from meshQualPlot import meshQualPlot
from funs import getMeshQualParams,getPlotData

 
# Setting directory to find .su2 files
os.chdir('c:\\Users\\floyd\\git\\Mesh-Deformation-RBF-Interpolation\\Classical RBF\\Meshes')

fNameInit = '/LS89_turbine.su2'

[f_init,v_init,elemType,bdryPnts_init] = getPlotData(fNameInit)

#[f,v,elemType,_] = getPlotData('/25x25.su2')  
#[f_init, v_init, elemType, bdryPnts_init] = getPlotData('/mesh_NACA0012_inv.su2')
#[f_init,v_init,elemType,bdryPnts_init] = getPlotData('/turbine_row.su2')
#[f,v,elemType] = getPlotData('/su2mesh.su2')
#[f,v,elemType] = getPlotData('/mesh_original.su2')
[alphas_0,_,_,_] = getMeshQualParams(f_init,v_init,elemType)





# Provide filenames of the meshes
#fileNames = ['/25x25_def_ps_none_noCorrect_DataRed_e1.su2','/25x25_def_ps_none_noCorrect_DataRed_e2.su2','/25x25_def_ps_none_noCorrect_DataRed_e3.su2','/25x25_def_ps_none_noCorrect_DataRed_e4.su2','/25x25_def_ps_none_correct_DataRed_e1.su2','/25x25_def_ps_none_correct_DataRed_e2.su2','/25x25_def_ps_none_correct_DataRed_e3.su2','/25x25_def_ps_none_correct_DataRed_e4.su2']
#fileNames = ['/25x25_def.su2']
#fileNames = ['/mesh_NACA0012_inv_def.su2']
fileNames = ['/LS89_turbine_def.su2']
#fileNames= ['/25x25_def_ref.su2','/25x25_def_err1.su2','/25x25_def_err2.su2','/25x25_def_err3.su2','/25x25_def_err4.su2','/25x25_def_err5.su2']
#fileNames = ['/mesh_NACA0012_inv_def_none_ref.su2','/mesh_NACA0012_inv_def_none_greedy.su2','/mesh_NACA0012_inv_def_ps_ref.su2','/mesh_NACA0012_inv_def_ps_greedy.su2']
#fileNames = ['/su2mesh_def_ps_ref.su2']
#fileNames = ['/mesh_original.su2']
#fileNames = ['/turbine_row_def.su2']

#graphNames = ['non-periodic', 'periodic in y', 'periodic in y,\nmoving boundaries, fixed vertices', 'periodic in y,\nmoving boundaries, moving vertices']
#graphNames = ['','Regular RBF', 'Greedy 1e-1', 'Greedy 1e-2', 'Greedy 1e-3', 'Greedy 1e-4', 'Greedy 1e-5']
graphNames = ['','','','']


#%% Mesh Quality plots
#meshQualPlot(fileNames,fNameInit,graphNames,alphas_0)


#%% MAKING A SCATTER PLOT OF THE BOUNDARY POINTS

### with the initial mesh ###
bdryScatterFuns.bdryScatter2(v_init,bdryPnts_init, fileNames[0])

### without the initial mesh ###
#bdryScatterFuns.bdryScatter(fileNames[0])


#%%
#plt.close('all')
#t = [2.2425,0.0140,0.1804,0.5408,1.6593,3.9731]
#N = [112,1,12,25,48,74]
#tol = [0,1,2,3,4,5]
#plt.figure()
#plt.plot(tol,t, '-x')
#plt.ylabel('time [sec]')
#plt.xlabel('Error tolerance [1e-x]')
#
#plt.figure()
#plt.plot(tol,N, '-x')
#plt.ylabel('Number of control nodes [-]')
#plt.xlabel('Error tolerance [1e-x]')



#%%
#plt.close('all')
#init = [19.8423,2.50666]
#new = [18.7065, 11.4972]
#mid = [16.8532, 10.6954]
#midx = [19.607,18.9836,18.0608,16.8532,15.3799,13.6639 ]
#midy =[3.74023,  6.16815,   8.49878,   10.6954,   12.7233,  14.5506]
#segmentx = [19.8423,19.3717 ,18.5955,17.5261,16.1803,14.5794]
#segmenty = [ 2.50666,4.9738,7.36249,9.63507,11.7557,13.6909]
#
#vec = [0.844329, 0.535825]
#a = 1
#midVecx = [16.8532, 16.8532+vec[0]*a]
#midVecy = [10.6954, 10.6954+a*vec[1]]
#   
#delta = [ 0.169354, -0.266861]
#plt.figure()
#plt.scatter(init[0],init[1])
#plt.scatter(new[0],new[1])
#plt.scatter(17.0226, 10.4285)
#plt.scatter(midx,midy)
#plt.scatter(mid[0],mid[1])
#plt.plot(midVecx,midVecy)
#plt.plot(segmentx,segmenty)
#ax = plt.gca()
#ax.set_aspect('equal')
#plt.show()
#
#print(new[0]-mid[0])
#print(new[1]-mid[1])