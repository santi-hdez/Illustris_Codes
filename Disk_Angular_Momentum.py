#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 17 17:31:00 2022

@author: santiago
"""

"""
                                                                              
 This procedure calculates the galaxy disk specific angular momentum by considering a dynamical decomposition in which the
 z-components of the specific angular momentums of star particles co-rotating and on nearly circular orbits are selected
 to compute the total angular momentum of the thin/cold disk. Only stars inside 2*half-mass stellar radius are considered.
 
 In addition, a representation of the z-components of the specific angular momentums of star particles shown 
 as a function of the specific binding energy is generated (orange dots). Theoretical angular momentums of stars 
 on perfect circular orbits for each binding energy is also represented (blue dots). Binding energies are
 computed relative to the total mass of the system within the virial radius, R200. 
 
 
 Execution: Jdisk = disk_angular_momentum(snap, ID, Mvirial200, simul, circ_limit)

 Input:

  simul:         Simulation name. 
                 Options: TNG300-1 --- TNG300-1 simulation.
                          TNG100-1 --- TNG100-1 simulation.
                          TNG50-1 --- TNG50-1 simulation.   
  snap:          Snapshot number in illustris simulation. It determines the redshift.      
  ID:            Subhalo ID in illustris simulation. 
  Mvirial200:    Total mass of the system within the virial radius, R200. (The corresponding halo ID in Illustris
                 must be searched. Field: 'Group_M_Crit200'). 
  circ_limit:    Circularity = Jz/Jcirc. Lower limit in the circularity for the selection of stars in the computation
                 of the disk angular momentum.                                                   
    

 Output:    
                                            
  Jdisk:       Angular momentum of the thin/cold galactic disk.
  
    
 * All quantities in Illustris units.                                                      

"""


def disk_angular_momentum(snap, ID, Mvirial200, simul='TNG100-1', circ_limit=0.8):
    
    import matplotlib.pyplot as plt
    import h5py
    import numpy as np
    import sys

    sys.path.append('/home/santiago/Documentos/SHPython/SHillustris')
                            
    from Illustris_API import get
    from Galaxy_Rotation import galaxy_rotation
    
    snap=str(snap)
    ID=str(ID)

    galaxy=get('http://www.illustris-project.org/api/'+simul+'/snapshots/'+snap+'/subhalos/'+ID)

    params_dam={'stars':'Coordinates,Velocities,GFM_StellarPhotometrics,Potential','gas':'Coordinates,Velocities,Masses,StarFormationRate,InternalEnergy,ElectronAbundance','dm':'Coordinates,Velocities'}

    cutout_dam=get('http://www.illustris-project.org/api/'+simul+'/snapshots/'+snap+'/subhalos/'+ID + "/cutout.hdf5",params_dam)

    
    with h5py.File(cutout_dam) as f_dam:
        
        G=43008.54006  #Universal gravitational constant (Illustris units)
        
        Jstars=[]
        
        Estars=[]
        
        Jcirc=[]
                
        J1=np.empty((0,3))
        
        CoordRotstars, CoordRotgas, VelRotstars, VelRotgas = galaxy_rotation(snap, ID, simul)
                
        for i in range(len(f_dam['PartType4']['Coordinates'][:,0])):
            
            #Radius of the rotated position of the star (see Galaxy_Rotation module):
            
            Ri = np.linalg.norm(CoordRotstars[:,i])
            
            #Mass around which the star orbits:
            
            M_orbitstars_i=-(f_dam['PartType4']['Potential'][i]*Ri)/G
            
            #Major semiaxis of the orbit:
            
            a_stars_i=G*M_orbitstars_i*Ri/(2*G*M_orbitstars_i-Ri*np.linalg.norm(VelRotstars[:,i])**2)
            
            #Binding energy of the star:
            
            Ebind200i=-G*Mvirial200/(a_stars_i*2)  #Potential + cinetic energy. 
            
            Estars=np.append(Estars,[Ebind200i],axis=0)
            
            #Angular momentum of the star:
            
            Ji=np.cross(CoordRotstars[:,i], VelRotstars[:,i])
            
            Jstars=np.append(Jstars,[Ji[2]],axis=0)
            
            #Theoretical angular momentum if the orbit were circular:
            
            Jcirc_i=np.sqrt(G*M_orbitstars_i*a_stars_i)
            
            Jcirc=np.append(Jcirc,[Jcirc_i],axis=0)
            
            #Calculation of the disk angular momentum:
            
            if Ri < 2*galaxy['halfmassrad_stars']:
            
                circ=Ji[2]/Jcirc_i
    
                if circ > circ_limit:
    
                    J1=np.append(J1,[Ji],axis=0)

    JX=sum(J1[:,0])
        
    JY=sum(J1[:,1])
        
    JZ=sum(J1[:,2])
        
    Jdisk=np.array([JX,JY,JZ])
    
    plt.scatter(Estars,Jcirc, s=1, marker='o')
    plt.scatter(Estars,Jstars, s=1, marker='v')
    
    plt.tick_params(labelsize=15)
    plt.xlabel('$E$ [km^{2}/s^{2}]', fontsize=15)
    plt.ylabel('$J$ [ckpc*km/h*s]', fontsize=15)
    
    
    return Jdisk, plt.show()

