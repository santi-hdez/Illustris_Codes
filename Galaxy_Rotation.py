#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 16:35:20 2022

@author: santiago
"""

"""
                                                                              
 This procedure rotates star and gas particles coordinates and velocities of the specified subhalo such that 
 if axes=[0,1] we have the face-on profile of the galaxy. For gas coordinates, an upper limit for the gas
 temperature can be imposed.

 Execution: CoordRotstars, CoordRotgas, VelRotstars, VelRotgas = galaxy_rotation(snap, ID, simul, T)

 Input:

  simul:       Simulation name. 
               Options: TNG300-1 --- TNG300-1 simulation.
                        TNG100-1 --- TNG100-1 simulation.
                        TNG50-1 --- TNG50-1 simulation.
  snap:        Snapshot number in illustris simulation. It determines the redshift.      
  ID:          Subhalo ID in illustris simulation.    
  T:           Upper limit for the gas temperature in kelvin.                                               
    

 Output:    
                                            
  CoordRotstars:   Rotated stars coordinates.
  CoordRotgas:     Rotated gas coordinates.
  VelRotstars:     Rotated stars velocities.
  VelRotgas:       Rotated gas velocities.

 * All quantities in Illustris units.                                                                 

"""


def galaxy_rotation(snap, ID, simul='TNG100-1', T='None'):

    import h5py
    import numpy as np
    import sys

    sys.path.append('/home/santiago/Documentos/SHPython/SHillustris')
                        
    from Illustris_API import get
    from Disk_Inertia_Tensor import disk_inertia_tensor
    from Cold_Gas import cold_gas

    snap=str(snap)
    ID=str(ID)

    galaxy=get('http://www.illustris-project.org/api/'+simul+'/snapshots/'+snap+'/subhalos/'+ID)

    params_rot={'stars':'Coordinates,Velocities,GFM_StellarPhotometrics,Potential','gas':'Coordinates,Velocities,Masses,StarFormationRate,InternalEnergy,ElectronAbundance','dm':'Coordinates,Velocities'}

    cutout_rot=get('http://www.illustris-project.org/api/'+simul+'/snapshots/'+snap+'/subhalos/'+ID + "/cutout.hdf5",params_rot)

    I, Mrot = disk_inertia_tensor(snap, ID, simul)

    with h5py.File(cutout_rot) as f_rot:
        
        CoordRotstars=[[],[],[]]
        
        CoordRotgas=[[],[],[]]
        
        VelRotstars=[[],[],[]]
        
        VelRotgas=[[],[],[]]
        
        GalPos=np.array([galaxy['pos_x'],galaxy['pos_y'],galaxy['pos_z']])
        
        
        #Rotation of coordinates:
            
        for i in range(len(f_rot['PartType4']['Coordinates'][:,0])):
            
            starparticle=np.dot(Mrot, np.array([[f_rot['PartType4']['Coordinates'][i,0]],[f_rot['PartType4']['Coordinates'][i,1]],[f_rot['PartType4']['Coordinates'][i,2]]])-np.transpose([GalPos]))
            
            CoordRotstars=np.append(CoordRotstars, starparticle, axis=1)
            
        
        if T == 'None':
            
            for i in range(len(f_rot['PartType0']['Coordinates'][:,0])):
            
                gasparticle=np.dot(Mrot, np.array([[f_rot['PartType0']['Coordinates'][i,0]],[f_rot['PartType0']['Coordinates'][i,1]],[f_rot['PartType0']['Coordinates'][i,2]]])-np.transpose([GalPos]))
            
                CoordRotgas=np.append(CoordRotgas, gasparticle, axis=1)
            
        else:
            
            ColdGasCoordinates, ColdGasMasses = cold_gas(snap, ID, T, simul)
            
            for i in range(len(ColdGasCoordinates[:,0])):
                
                gasparticle=np.dot(Mrot, np.array([[ColdGasCoordinates[i,0]],[ColdGasCoordinates[i,1]],[ColdGasCoordinates[i,2]]])-np.transpose([GalPos]))
                
                CoordRotgas=np.append(CoordRotgas, gasparticle, axis=1)
            
    
        #Rotation of velocities:
            
        for i in range(len(f_rot['PartType4']['Velocities'][:,0])):
        
            starparticle=np.dot(Mrot, np.array([[f_rot['PartType4']['Velocities'][i,0]],[f_rot['PartType4']['Velocities'][i,1]],[f_rot['PartType4']['Velocities'][i,2]]]))
        
            VelRotstars=np.append(VelRotstars, starparticle, axis=1)
        
        for i in range(len(f_rot['PartType0']['Velocities'][:,0])):
        
            gasparticle=np.dot(Mrot, np.array([[f_rot['PartType0']['Velocities'][i,0]],[f_rot['PartType0']['Velocities'][i,1]],[f_rot['PartType0']['Velocities'][i,2]]]))
        
            VelRotgas=np.append(VelRotgas, gasparticle, axis=1)
            
            
    return CoordRotstars, CoordRotgas, VelRotstars, VelRotgas
        
        
        
        
        
        
        
        
        
        
        