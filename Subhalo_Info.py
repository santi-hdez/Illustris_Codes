#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 20:11:18 2022

@author: santiago
"""

"""
                                                                              
 This procedure gives the general information about a specific subhalo (mass, position, number of particles, components, etc). 

 Execution: info = subhalo_info(snap, ID, simul)

 Input:

  simul:       Simulation name. 
               Options: TNG300-1 --- TNG300-1 simulation.
                        TNG100-1 --- TNG100-1 simulation.
                        TNG50-1 --- TNG50-1 simulation.
  snap:        Snapshot number in illustris simulation. It determines the redshift.      
  ID:          Subhalo ID in illustris simulation.                                                   
    

 Output:    
                                            
  info:        General information about the specified subhalo.
  
  * All quantities in Illustris units.
                                                                  

"""


def subhalo_info(snap, ID, simul='TNG100-1'):
    
    import sys

    sys.path.append('/home/santiago/Documentos/SHPython/SHillustris')
                            
    from Illustris_API import get
    
    snap=str(snap)
    ID=str(ID)

    info=get('http://www.illustris-project.org/api/'+simul+'/snapshots/'+snap+'/subhalos/'+ID)

    return info












