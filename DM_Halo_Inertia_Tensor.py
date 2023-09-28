#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 16:38:13 2022

@author: santiago
"""

"""
                                                                              
 This procedure computes the dark matter halo inertia tensor considering all dark matter particles of the subhalo.
 It also provides a matrix with the dark matter halo principal axes of inertia.

 Execution: IDm, IDm_main_axes = dm_halo_inertia_tensor(snap, ID, simul)

 Input:

  simul:       Simulation name. 
               Options: TNG300-1 --- TNG300-1 simulation.
                        TNG100-1 --- TNG100-1 simulation.
                        TNG50-1 --- TNG50-1 simulation.
  snap:        Snapshot number in illustris simulation. It determines the redshift.      
  ID:          Subhalo ID in illustris simulation.                                                   
    

 Output:    
                                            
  IDm:             Inertia tensor of the dark matter halo.
  Dm_main_axes:    Matrix with the dark matter halo principal axes of inertia (eigenvectors).
                                                                  
 * All quantities in Illustris units.

"""

def dm_halo_inertia_tensor(snap, ID, simul='TNG100-1'):

    import h5py
    import numpy as np
    import sys

    sys.path.append('/home/santiago/Documentos/SHPython/SHillustris')
                        
    from Illustris_API import get

    snap=str(snap)
    ID=str(ID)


    galaxy=get('http://www.illustris-project.org/api/'+simul+'/snapshots/'+snap+'/subhalos/'+ID)

    params_dmit={'stars':'Coordinates,Velocities,GFM_StellarPhotometrics,Potential','gas':'Coordinates,Velocities,Masses,StarFormationRate,InternalEnergy,ElectronAbundance','dm':'Coordinates,Velocities'}

    cutout_dmit=get('http://www.illustris-project.org/api/'+simul+'/snapshots/'+snap+'/subhalos/'+ID + "/cutout.hdf5",params_dmit)


    GalPos=np.array([galaxy['pos_x'],galaxy['pos_y'],galaxy['pos_z']]) 

    Mdm=[]

    if simul == 'TNG300-1':

        Mdm_i =	0.00398342749867548

    if simul == 'TNG100-1':

        Mdm_i = 0.000505574296436975
    
    if simul == 'TNG50-1':

        Mdm_i = 0.0000307367708626464


    with h5py.File(cutout_dmit) as f_dmit:
        
    
        for i in range(len(f_dmit['PartType1']['Coordinates'][:,1])):
        
            Mdm=np.append(Mdm, [Mdm_i], axis=0)
        
        
        CoordDm=f_dmit['PartType1']['Coordinates'][:]
    

    for i in range(3):
        
        CoordDm[:,i]-=GalPos[i]
        

    IDm=np.zeros((3,3), dtype='float32')

    IDm[0,0]=np.sum(Mdm*(CoordDm[:,1]*CoordDm[:,1]+CoordDm[:,2]*CoordDm[:,2]))
    IDm[1,1]=np.sum(Mdm*(CoordDm[:,0]*CoordDm[:,0]+CoordDm[:,2]*CoordDm[:,2]))
    IDm[2,2]=np.sum(Mdm*(CoordDm[:,0]*CoordDm[:,0]+CoordDm[:,1]*CoordDm[:,1]))
    IDm[0,1]=-1*np.sum(Mdm*(CoordDm[:,0]*CoordDm[:,1]))
    IDm[0,2]=-1*np.sum(Mdm*(CoordDm[:,0]*CoordDm[:,2]))
    IDm[1,2]=-1*np.sum(Mdm*(CoordDm[:,1]*CoordDm[:,2]))
    IDm[1,0]=IDm[0,1]
    IDm[2,0]=IDm[0,2]
    IDm[2,1]=IDm[1,2]


    #Get eigenvalues and normalized right eigenvectors:

    eigen_valuesDm, Dm_main_axes = np.linalg.eig(IDm)

    #Sort ascending the eigenvalues:

    sort_index = np.argsort(eigen_valuesDm)
    eigen_valuesDm = eigen_valuesDm[sort_index]


    #Permute the eigenvectors into this order, which is the rotation matrix which orients the
    #principal axes to the cartesian x,y,z axes, such that if axes=[0,1] we have the face-on profile:

    Dm_main_axes = np.matrix((Dm_main_axes[:,sort_index[0]],
                               Dm_main_axes[:,sort_index[1]],
                               Dm_main_axes[:,sort_index[2]]))


    Dm_main_axes=np.array(Dm_main_axes)
    
    return IDm, Dm_main_axes




