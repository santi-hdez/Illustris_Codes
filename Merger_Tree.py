#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 23:48:21 2022

@author: santiago
"""

"""
                                                                              
 This procedure gives the snapshots numbers and the subhalo IDs of the main progenitor branch. It can be selected between
 the sublink algorithm and the lhalotree algorithm.
 
 Execution: snaps, SubhaloIDs = merger_tree(ID, simul, algorithm, flag)

 Input:

  simul:       Simulation name. 
               Options: TNG300-1 --- TNG300-1 simulation.
                        TNG100-1 --- TNG100-1 simulation.
                        TNG50-1 --- TNG50-1 simulation.    
  ID:          Subhalo ID in illustris simulation at z=0.   
  algorithm:   Indicates which algortihm to use to obtain the main progenitor branch.
               Options: sublink --- sublink algorithm.
                        lhalotree --- lhalotree algorithm.
  flag:        Indicates if the output must contain the complete main progenitor branch or only those 'full'
               snapshots in illustris.
               Options: all --- all main progenitor branch.
                        full --- only 'full' snapshots.
               
                                                
 Output:    
                                            
  snaps:        Snapshot numbers of the main progenitor branch.
  SubhaloIDs:   Subhalo IDs of the main progenitor branch.
                                                                  

"""

def merger_tree(ID, simul='TNG100-1', algorithm='sublink', flag='all'):

    import h5py
    import numpy as np
    import sys

    sys.path.append('/home/santiago/Documentos/SHPython/SHillustris')
                            
    from Illustris_API import get
    
    np.set_printoptions(threshold=sys.maxsize)
    
    ID=str(ID)

    mpb=get('http://www.illustris-project.org/api/'+simul+'/snapshots/99/subhalos/'+ID+'/'+algorithm+'/mpb.hdf5')
        
    f = h5py.File(mpb)
     
    if algorithm == 'sublink':
        
        if flag == 'all':
            
            snaps = f['SnapNum'][:]
            
            SubhaloIDs = (f['SubfindID'][:])
            
        if flag == 'full':
            
            limit=len(f['SnapNum'][:])-1
            
            i=np.array([0,8,15,21,27,32,40,49,59,66,74,78,82,86,88,91,93,95,96,97])
            
            j=np.where(i<=limit)
            
            snaps = f['SnapNum'][i[j]]
            
            SubhaloIDs = (f['SubfindID'][i[j]])
            
    if algorithm == 'lhalotree':
           
        if flag == 'all':
               
            snaps = f['SnapNum'][:]
               
            SubhaloIDs = (f['SubhaloNumber'][:])
               
        if flag == 'full':
            
            limit=len(f['SnapNum'][:])-1
               
            i=np.array([0,8,15,21,27,32,40,49,59,66,74,78,82,86,88,91,93,95,96,97])
            
            j=np.where(i<=limit)
               
            snaps = f['SnapNum'][i[j]]
            
            SubhaloIDs = (f['SubhaloNumber'][i[j]])
        
        
    return snaps, SubhaloIDs

    
   
    