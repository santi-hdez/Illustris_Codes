#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 17:56:18 2022

@author: santiago
"""

"""
                                                                              
 This procedure creates a 3D interactive graphic of an specific subhalo in illustris simulation. *Valid in
 Spyder environment previously running the command '%matplotlib qt5' in the console. After the execution run
 '%matplotlib inline' to desactivate the interactive mode.
 
 Execution: interactive_galaxy_3d(snap, ID, simul, T, particles)

 Input:

  simul:       Simulation name. 
               Options: TNG300-1 --- TNG300-1 simulation.
                        TNG100-1 --- TNG100-1 simulation.
                        TNG50-1 --- TNG50-1 simulation.   
  snap:        Snapshot number in illustris simulation. It determines the redshift.      
  ID:          Subhalo ID in illustris simulation.    
  T:           Upper limit for the gas temperature in kelvin.     
  particles:   Particles to draw.
               Options: gas --- only represents gas particles.
                        stars --- only represents stars particles.
                        both --- represents both gas and stars particles.                                   
    

 Output:    
                                            
  A 3D interactive graphic of the subhalo.    

 * All quantities in Illustris units.                                                      

"""


def interactive_galaxy_3d(snap, ID, simul='TNG100-1', T='None', particles='gas'):
    
    import h5py
    import matplotlib.pyplot as plt
    import sys

    sys.path.append('/home/santiago/Documentos/SHPython/SHillustris')
                        
    from Illustris_API import get
    from Cold_Gas import cold_gas

    snap=str(snap)
    ID=str(ID)


    params_int={'stars':'Coordinates,Velocities,GFM_StellarPhotometrics,Potential','gas':'Coordinates,Velocities,Masses,StarFormationRate,InternalEnergy,ElectronAbundance','dm':'Coordinates,Velocities'}

    cutout_int=get('http://www.illustris-project.org/api/'+simul+'/snapshots/'+snap+'/subhalos/'+ID + "/cutout.hdf5",params_int)


    fig=plt.figure(figsize=(10,10))
    ax=fig.add_subplot(111,projection='3d')

    with h5py.File(cutout_int) as f_int:
        
        if T == 'None':
            
            xg=f_int['PartType0']['Coordinates'][:,0]
            yg=f_int['PartType0']['Coordinates'][:,1]
            zg=f_int['PartType0']['Coordinates'][:,2]
        
        else:
        
            ColdGasCoordinates, ColdGasMasses = cold_gas(snap, ID, T, simul)
            
            xg=ColdGasCoordinates[:,0]
            yg=ColdGasCoordinates[:,1]
            zg=ColdGasCoordinates[:,2] 
            
        xs=f_int['PartType4']['Coordinates'][:,0]
        ys=f_int['PartType4']['Coordinates'][:,1]
        zs=f_int['PartType4']['Coordinates'][:,2]
        
    if particles == 'gas':
        
        ax.scatter(xg, yg, zg, s=2, marker='.')
        
    if particles == 'stars':
        
        ax.scatter(xs, ys, zs, s=2, marker='.')
        
    if particles == 'both':     
        
        ax.scatter(xg, yg, zg, s=2, marker='.')
        ax.scatter(xs, ys, zs, s=2, marker='.')
        
    
    return plt.show()








