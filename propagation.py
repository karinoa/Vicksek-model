# -*- coding: utf-8 -*-
"""
Created on Wed May  2 15:47:28 2018

@author: Koen
"""
import numpy as np
import math
from core import *


def Eulerpropstep():#positions,velocities, more?):
    #update position: position += timestep*velocity
    #check nearest neighbours and check if boundary (defined by angle > pi)
    #calculate average angle
    #calculate force: F = Fself + Fboundary + Frepulsion
    #update velocity: velocity += timestep*force
    
    return #positions,velocities,angle?,forces
    
    
def Neighbourcheck():
    #check nearest neighbours
    #particles are neighbours if distance between centers is <= 2.7*MEAN_RADIUS
    return
    
def Boundarycheck():
    #check if boundary (defined by angle > pi)
    return
    
def Distancematrix(system):
    distances = np.zeroes(N_PARTICLES,N_PARTICLES)
    for i in range(N_PARTICLES): #for every particle
        for j in range(N_PARTICLES): #in relation to another particle
            if j != i:  #if they're not the same
                distances(i,j)= math.sqrt((system[i,COLUMN_REVERSE_MAPPING['x']]-system[j,COLUMN_REVERSE_MAPPING['x']])^2 + (system[i,COLUMN_REVERSE_MAPPING['y']]-system[j,COLUMN_REVERSE_MAPPING['y']])^2)
    print(distances)
                
    return distances