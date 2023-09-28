#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 21:55:22 2022

@author: santiago
"""

"""
                                                                              
 This procedure creates an isophotes representation of an specific subhalo in illustris simulation. It represents 
 the contours of stellar magnitudes in one of the eight available bands in a face-on profile of the subhalo.
 
 Execution: isophotes(snap, ID, simul, band, bins_, levels_, range_)

 Input:

  simul:       Simulation name. 
               Options: TNG300-1 --- TNG300-1 simulation.
                        TNG100-1 --- TNG100-1 simulation.
                        TNG50-1 --- TNG50-1 simulation.   
  snap:        Snapshot number in illustris simulation. It determines the redshift.      
  ID:          Subhalo ID in illustris simulation.   
  band:        Stellar magnitudes to represent. Eight bands. 
               Options: U --- 0
                        B --- 1
                        V --- 2
                        K --- 3
                        g --- 4
                        r --- 5
                        i --- 6
                        z --- 7
  bins_:       Bins for the isophotes representation. 
  levels_:     Contour levels for the isophotes representation.                                                 
  range_:      Range in x and y axis for the representation.  

 Output:    
                                            
 Isophotes representation of the subhalo. Face-on profile.

 * All quantities in Illustris units.                                                      

"""

def isophotes(snap, ID, simul='TNG100-1', band=2, bins_=[100,100], levels_=30, range_=[[-20,20],[-20,20]]):

    import matplotlib.pyplot as plt
    import h5py
    import sys

    sys.path.append('/home/santiago/Documentos/SHPython/SHillustris')
                        
    from Illustris_API import get
    from Galaxy_Rotation import galaxy_rotation
    
    snap=str(snap)
    ID=str(ID)

    params_isoph={'stars':'Coordinates,Velocities,GFM_StellarPhotometrics,Potential','gas':'Coordinates,Velocities,Masses,StarFormationRate,InternalEnergy,ElectronAbundance','dm':'Coordinates,Velocities'}

    cutout_isoph=get('http://www.illustris-project.org/api/'+simul+'/snapshots/'+snap+'/subhalos/'+ID + "/cutout.hdf5",params_isoph)
   
    CoordRotstars, CoordRotgas, VelRotstars, VelRotgas = galaxy_rotation(snap, ID, simul)

    with h5py.File(cutout_isoph) as f_isoph:
        
        x=CoordRotstars[0,:] 
        y=CoordRotstars[1,:] 
        dens=10**(-0.4*f_isoph['PartType4']['GFM_StellarPhotometrics'][:,band])  
    
    fig = plt.figure(figsize=(8,8))
    
    plt.subplot(121)
    
    counts,xbins,ybins,image=plt.hist2d(x, y, weights=dens, bins=bins_, cmap='gist_heat_r', range=range_)
    
    plt.subplot(111)
    
    plt.contour(counts.transpose(), extent=[xbins[0],xbins[-1],ybins[0],ybins[-1]], linewidths=1,
                levels=levels_, cmap='rainbow')

    plt.tick_params(labelsize=15)
    
    plt.xlabel('$\Delta x$ [ckpc/h]', fontsize=15)
    
    plt.ylabel('$\Delta y$ [ckpc/h]', fontsize=15)
    
    return plt.show()




