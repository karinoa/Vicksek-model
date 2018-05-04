# -*- coding: utf-8 -*-
"""
Created on Wed May  2 15:47:28 2018

@author: Koen
"""
import numpy as np
from core import *



def Eulerpropstep():#positions,velocities, more?):
    #update position: position += timestep*velocity
    #check nearest neighbours and check if boundary (defined by angle > pi)
    #calculate average angle
    #calculate force: F = Fself + Fboundary + Frepulsion
    #update velocity: velocity += timestep*force
    
    return #positions,velocities,angle?,forces
    
    
def get_neighbours(system): #check nearest neighbours
    neighbours = np.zeros(shape=(N_PARTICLES,N_PARTICLES), dtype = bool) #particles are neighbours if distance between centers is <= 2.7*MEAN_RADIUS
    distances = get_distances(system)    
    for i in range(N_PARTICLES):
        neighbours[i,i] = 'TRUE'
        for j in range(i,N_PARTICLES):
            if distances[i,j]<=2.7*MEAN_RADIUS:
                neighbours[i,j] = 'TRUE'
    return neighbours
    
def Boundarycheck():
    #check if boundary (defined by angle > pi)
    return
    
def get_distances(system):
    distancematrix = np.zeros(shape=(N_PARTICLES,N_PARTICLES))
    for i in range(N_PARTICLES): #for every particle
        distancematrix[i,i] = 0     #by definition
        for j in range(N_PARTICLES): #in relation to another particle
            if j != i:  #if they're not the same    
                distancematrix[i,j]= np.sqrt((system[i,COLUMN_REVERSE_MAPPING['x']]-system[j,COLUMN_REVERSE_MAPPING['x']])**2 + (system[i,COLUMN_REVERSE_MAPPING['y']]-system[j,COLUMN_REVERSE_MAPPING['y']])**2)
    return distancematrix
    
print(get_neighbours(system)[0])