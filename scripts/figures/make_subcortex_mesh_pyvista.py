#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 21:37:52 2020

@author: sidchopra
"""


import numpy as np
import pyvista as pv
from numpy import inf


#import matplotlib.pyplot as plt


#max_degree = np.array([0,63, 63, 58, 58, 41, 41, 6, 6, 21, 21, 12, 12, 5, 5])
max_degree = np.array([0,63, 63, 58, 58,21, 21, 6, 6, 21, 21, 12, 12, 5, 5]) # medication effect max halfed to stop influence of outlier
for i in range(1,15):
    mesh = pv.read('/Users/sidchopra/Dropbox/Sid/R_files/STAGES_fmri/output/subcortex/tian_lh.vtk')
    meshSmooth = mesh.smooth(n_iter=200)
    
    
    ##get scalars
    cells = mesh.active_scalars
    unique_elements, counts_elements = np.unique(cells, return_counts=True)
    
    #read in degree file
    degreefile ='/Users/sidchopra/Dropbox/Sid/R_files/STAGES_fmri/output/subcortex/sub_degree' + str(i) + '.txt'
    degree = np.genfromtxt(degreefile, delimiter=',')
    lh_degree = degree[0:16]
    scalarsDegree = np.repeat(lh_degree, counts_elements, axis=0)
    
    
    meshSmooth.plot(scalars=scalarsDegree, cmap="Reds", 
          background="White", 
          parallel_projection=True, 
          clim = [0,max_degree[i]],
          cpos=[1.5, 0.5, -1], screenshot='/Users/sidchopra/Dropbox/Sid/R_files/STAGES_fmri/output/subcortex/imgs/sub_'+str(i)+'_l.png')
    
    
    mesh2 = pv.read('/Users/sidchopra/Dropbox/Sid/R_files/STAGES_fmri/output/subcortex/tian_rh.vtk')
    mesh2Smooth = mesh2.smooth(n_iter=200)
    rh_degree = degree[16:32]
    scalarsDegree2 = np.repeat(rh_degree, counts_elements, axis=0)
    
    
    mesh2Smooth.plot(scalars=scalarsDegree2, cmap="Reds", 
                     background="White", 
                     parallel_projection=True, 
                     clim = [0,max_degree[i]],
                     cpos=[-10, 4, -6],
                     show_scalar_bar=True, screenshot='/Users/sidchopra/Dropbox/Sid/R_files/STAGES_fmri/output/subcortex/imgs/sub_'+str(i)+'_r.png')
    
    
    
    
    






# Make legend 
mesh = pv.read('/Users/sidchopra/Dropbox/Sid/R_files/STAGES_fmri/output/subcortex/tian_lh.vtk')
meshSmooth = mesh.smooth(n_iter=200)


##get scalars
cells = mesh.active_scalars
unique_elements, counts_elements = np.unique(cells, return_counts=True)

    #read in degree file
degreefile ='/Users/sidchopra/Dropbox/Sid/R_files/STAGES_fmri/output/subcortex/sub_degree' + str(i) + '.txt'
degree = np.genfromtxt(degreefile, delimiter=',')
lh_degree = degree[0:16]
scalarsDegree = np.repeat(lh_degree, counts_elements, axis=0)


meshSmooth.plot(scalars=scalarsDegree, cmap="Reds", 
                background="White", 
                parallel_projection=True, 
                clim = [0,max_degree[i]],
                cpos=[1.5, 0.5, -1], screenshot='/Users/sidchopra/Dropbox/Sid/R_files/STAGES_fmri/output/subcortex/imgs/sub_'+str(i)+'_l.png')


mesh2 = pv.read('/Users/sidchopra/Dropbox/Sid/R_files/STAGES_fmri/output/subcortex/tian_rh.vtk')
mesh2Smooth = mesh2.smooth(n_iter=200)
rh_degree = degree[16:32]
scalarsDegree2 = np.repeat(rh_degree, counts_elements, axis=0)


mesh2Smooth.plot(scalars=scalarsDegree2, cmap="Reds", 
                     background="White", 
                     parallel_projection=True, 
                     clim = [0,max_degree[i]],
                     cpos=[-10, 4, -6],
                     show_scalar_bar=True, screenshot='/Users/sidchopra/Dropbox/Sid/R_files/STAGES_fmri/output/subcortex/imgs/sub_'+str(i)+'_r.png')

