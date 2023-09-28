#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 18:22:00 2022

@author: santiago
"""

"""
                                                                              
 This procedure selects the coordinates and masses of the galaxy cold gas under a certain temperature. 
 
 Execution: ColdGasCoordinates, ColdGasMasses = cold_gas(snap, ID, T, simul)

 Input:

  simul:       Simulation name. 
               Options: TNG300-1 --- TNG300-1 simulation.
                        TNG100-1 --- TNG100-1 simulation.
                        TNG50-1 --- TNG50-1 simulation.  
  snap:        Snapshot number in illustris simulation. It determines the redshift.      
  ID:          Subhalo ID in illustris simulation.   
  T:           Upper limit for the gas temperature in kelvin.                                                 
    

 Output:    
                                            
  ColdGasCoordinates:   Coordinates of gas particles with temperature under T.
  ColdGasMasses:        Masses of gas particles with temperature under T.    

 * All quantities in Illustris units.                                                             

"""


def cold_gas(snap, ID, T, simul='TNG100-1'):
    
    import h5py
    import numpy as np
    import sys     

    sys.path.append('/home/santiago/Documentos/SHPython/SHillustris')  

    from Illustris_API import get                
     
    snap=str(snap)
    ID=str(ID)

    params_cold={'stars':'Coordinates,Velocities,GFM_StellarPhotometrics,Potential','gas':'Coordinates,Velocities,Masses,StarFormationRate,InternalEnergy,ElectronAbundance','dm':'Coordinates,Velocities'}

    cutout_cold=get('http://www.illustris-project.org/api/'+simul+'/snapshots/'+snap+'/subhalos/'+ID + "/cutout.hdf5",params_cold)

    with h5py.File(cutout_cold) as f_cold:
    

        k=1.380650*(10**(-16)) #(CGS)
        mp=1.6726216*(10**(-24)) #(CGS)
        Xh=0.76
    
        ColdGasCoordinates=np.empty((0,3))
    
        ColdGasMasses=[]
    
        for i in range(len(f_cold['PartType0']['Coordinates'][:,0])):
        
            mu=(4*mp)/(1+3*Xh+4*Xh*f_cold['PartType0']['ElectronAbundance'][i])
    
            Ti=(2*f_cold['PartType0']['InternalEnergy'][i]*mu*10**10)/(3*k) 
          
            if Ti<T:
            
                ColdGasCoordinates=np.append(ColdGasCoordinates,[f_cold['PartType0']['Coordinates'][i,:]], axis=0)
            
                ColdGasMasses=np.append(ColdGasMasses,[f_cold['PartType0']['Masses'][i]], axis=0)
            
    return ColdGasCoordinates, ColdGasMasses
    









