#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 02:48:06 2022

@author: santiago
"""

"""
                                                                              
 This procedure computes the specific angular momentum of the dark matter halo considering all dark matter
 particles of the subhalo. 

 Execution: Jdm = dm_halo_angular_momentum(snap, ID, simul)

 Input:

  simul:       Simulation name. 
               Options: TNG300-1 --- TNG300-1 simulation.
                        TNG100-1 --- TNG100-1 simulation.
                        TNG50-1 --- TNG50-1 simulation.
  snap:        Snapshot number in illustris simulation. It determines the redshift.      
  ID:          Subhalo ID in illustris simulation.                                                   
    

 Output:    
                                            
  Jdm:        Dark matter halo angular momentum of the subhalo.
                           
 * All quantities in Illustris units.                                       

"""


def dm_halo_angular_momentum(snap, ID, simul='TNG100-1'):

    import h5py
    import numpy as np
    import sys

    sys.path.append('/home/santiago/Documentos/SHPython/SHillustris')
                            
    from Illustris_API import get
    
    snap=str(snap)
    ID=str(ID)

    params_dmam={'stars':'Coordinates,Velocities,GFM_StellarPhotometrics,Potential','gas':'Coordinates,Velocities,Masses,StarFormationRate,InternalEnergy,ElectronAbundance','dm':'Coordinates,Velocities'}

    cutout_dmam=get('http://www.illustris-project.org/api/'+simul+'/snapshots/99/subhalos/'+ID + "/cutout.hdf5",params_dmam)


    with h5py.File(cutout_dmam) as f_dmam:
    
        #Compute the angular momentum for each dm particle:
    
        Jtemp=np.empty((0,3))
    
        for i in range(len(f_dmam['PartType1']['Coordinates'][:,0])):
        
            Ji=np.cross(f_dmam['PartType1']['Coordinates'][i,:], f_dmam['PartType1']['Velocities'][i,:])
        
            Jtemp=np.append(Jtemp,[Ji],axis=0)
        
        #Total angular momentum in each axis:
    
        JX=sum(Jtemp[:,0])
    
        JY=sum(Jtemp[:,1])
    
        JZ=sum(Jtemp[:,2])
    
        Jdm=np.array([JX,JY,JZ])
        
    return Jdm








