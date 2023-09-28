#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 16:16:05 2022

@author: santiago
"""

"""
                                                                              
 This procedure calculates the inertia tensor of a galaxy disk by considering the gas cells inside 2*half-mass stellar 
 radius and with active star formation rate. It also provides a rotation matrix for edge-on or face-on representations of
 the galaxy.
 
 Execution: I, Mrot = disk_inertia_tensor(snap, ID, simul)

 Input:

  simul:       Simulation name. 
               Options: TNG300-1 --- TNG300-1 simulation.
                        TNG100-1 --- TNG100-1 simulation.
                        TNG50-1 --- TNG50-1 simulation.   
  snap:        Snapshot number in illustris simulation. It determines the redshift.      
  ID:          Subhalo ID in illustris simulation.                                                    
    

 Output:    
                                            
  I:           Inertia tensor of the galaxy disk.
  Mrot:        Rotation matrix that contains the inertia tensor eigenvectors. It orients the 
               principal axes to the cartesian x,y,z axes, such that if axes=[0,1] we have the
               face-on profile of the galaxy.     

 * All quantities in Illustris units.                                                      

"""

def disk_inertia_tensor(snap, ID, simul='TNG100-1'):

    import h5py
    import numpy as np
    import sys

    sys.path.append('/home/santiago/Documentos/SHPython/SHillustris')
                            
    from Illustris_API import get
    
    snap=str(snap)
    ID=str(ID)

    galaxy=get('http://www.illustris-project.org/api/'+simul+'/snapshots/'+snap+'/subhalos/'+ID)

    params_dit={'stars':'Coordinates,Velocities,GFM_StellarPhotometrics,Potential','gas':'Coordinates,Velocities,Masses,StarFormationRate,InternalEnergy,ElectronAbundance','dm':'Coordinates,Velocities'}

    cutout_dit=get('http://www.illustris-project.org/api/'+simul+'/snapshots/'+snap+'/subhalos/'+ID + "/cutout.hdf5",params_dit)



    with h5py.File(cutout_dit) as f_dit:
    
        rHalf=galaxy['halfmassrad_stars']
    
        GalPos=np.array([galaxy['pos_x'],galaxy['pos_y'],galaxy['pos_z']]) 
    
        dist_gas=[]
    
        for i in range(len(f_dit['PartType0']['Coordinates'][:,1])):
    
            dist=np.linalg.norm(f_dit['PartType0']['Coordinates'][i,:]-GalPos)
        
            dist_gas=np.append(dist_gas, [dist], axis=0)
    
        Gas_index=np.where((dist_gas <= 2.0*rHalf) & (f_dit['PartType0']['StarFormationRate'][:]>0.0))
    
        Gas_index=np.squeeze(Gas_index)
    
        Masses=f_dit['PartType0']['Masses'][Gas_index]
        
        Coord=f_dit['PartType0']['Coordinates'][Gas_index,:]

        for i in range(3):
        
            Coord[:,i]-=GalPos[i]

        I=np.zeros((3,3), dtype='float32')

        I[0,0]=np.sum(Masses*(Coord[:,1]*Coord[:,1]+Coord[:,2]*Coord[:,2]))
        I[1,1]=np.sum(Masses*(Coord[:,0]*Coord[:,0]+Coord[:,2]*Coord[:,2]))
        I[2,2]=np.sum(Masses*(Coord[:,0]*Coord[:,0]+Coord[:,1]*Coord[:,1]))
        I[0,1]=-1*np.sum(Masses*(Coord[:,0]*Coord[:,1]))
        I[0,2]=-1*np.sum(Masses*(Coord[:,0]*Coord[:,2]))
        I[1,2]=-1*np.sum(Masses*(Coord[:,1]*Coord[:,2]))
        I[1,0]=I[0,1]
        I[2,0]=I[0,2]
        I[2,1]=I[1,2]
        

    #Get eigenvalues and normalized right eigenvectors:

    eigen_values, rotation_matrix = np.linalg.eig(I)

    #Sort ascending the eigenvalues:

    sort_index = np.argsort(eigen_values)
    eigen_values = eigen_values[sort_index]


    #Permute the eigenvectors into this order, which is the rotation matrix which orients the
    #principal axes to the cartesian x,y,z axes, such that if axes=[0,1] we have the face-on profile:

    new_matrix = np.matrix((rotation_matrix[:,sort_index[0]],
                            rotation_matrix[:,sort_index[1]],
                            rotation_matrix[:,sort_index[2]]))



    Mrot=np.array(new_matrix)
    
    return I, Mrot









