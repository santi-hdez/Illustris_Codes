#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 16:29:08 2022

@author: santiago
"""

"""
                                                                              
 This procedure allows to use the Illustris API given an API_key.                                                                                                                                              


  API_key:   Api key associated with an specific account in Illustris.
                                                     
  
"""


def get(path, params=None):
    
    import requests
    
    API_key='6e581b5220d62407184d08c36f9579a5'
    
    # make HTTP GET request to path
    r = requests.get(path, params=params, headers={"api-key":API_key})

    # raise exception if response code is not HTTP SUCCESS (200)
    r.raise_for_status()

    if r.headers['content-type'] == 'application/json':
        return r.json() # parse json responses automatically
    
    if 'content-disposition' in r.headers:
        filename = r.headers['content-disposition'].split("filename=")[1]
        with open(filename, 'wb') as f:
            f.write(r.content)
        return filename # return the filename string

    return r
    
    
    
    
    