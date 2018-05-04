# -*- coding: utf-8 -*-
"""
Created on Wed May  2 15:47:28 2018

@author: Koen
"""
import numpy as np
from core import *

K = 2
NEIGHBOUR_CUTOFF = 2.7 * MEAN_RADIUS
def Eulerpropstep(system):#positions,velocities, more?):
    distances,directions = get_distances(system)    #get distances between particles and the direction of that vector
    neighbours = get_neighbours(system,distances)   #get neighbours of each particle
    F = get_forces(system,distances, directions)
    #update position: position += timestep*velocity
    #check nearest neighbours and check if boundary (defined by angle > pi)
    #calculate average angle
    #calculate force: F = Fself + Fboundary + Frepulsion
    #update velocity: velocity += timestep*force
    print(F[3])
    return #positions,velocities,angle?,forces
    
def get_forces(system,distances,directions):
    forcematrix = np.zeros(shape=(N_PARTICLES,N_PARTICLES,2))
    for i in range(N_PARTICLES):
        for j in range(N_PARTICLES):
            rsum = system[i,COLUMN_REVERSE_MAPPING['r']] + system[j,COLUMN_REVERSE_MAPPING['r']] #r1 + r2
            dr = distances[i,j] #distance between particle centres
            if i!=j and dr <= rsum:
                force_repulsion = -K* (rsum/dr - 1)*directions[i,j]

            else: force_repulsion = 0
    forcematrix += force_repulsion
    return forcematrix
    
def get_neighbours(system,distances): #check nearest neighbours
    neighbours = np.zeros(shape=(N_PARTICLES,N_PARTICLES), dtype = bool) #particles are neighbours if distance between centers is <= 2.7*MEAN_RADIUS    
    for i in range(N_PARTICLES):
        neighbours[i,i] = 'TRUE'
        for j in range(N_PARTICLES):
            if distances[i,j] <= NEIGHBOUR_CUTOFF:
                neighbours[i,j] = 'TRUE'
    return neighbours
    
def Boundarycheck():
    #check if boundary (defined by angle > pi)
    return

def get_distances(system):
    distancematrix = np.zeros(shape=(N_PARTICLES,N_PARTICLES))
    directionmatrix = np.zeros(shape=(N_PARTICLES,N_PARTICLES,2))
    for i in range(N_PARTICLES): #for every particle
        distancematrix[i,i] = 0     #by definition
        directionmatrix[i,i] = 0
        for j in range(N_PARTICLES): #in relation to another particle
            if j != i:  #if they're not the same 
                dx = system[j,COLUMN_REVERSE_MAPPING['x']]-system[i,COLUMN_REVERSE_MAPPING['x']]
                dy = system[j,COLUMN_REVERSE_MAPPING['y']]-system[i,COLUMN_REVERSE_MAPPING['y']]
                distancematrix[i,j] = np.sqrt((dx)**2 + (dy)**2)
                dx = dx / distancematrix[i,j]
                dy = dy / distancematrix[i,j]
                directionmatrix[i,j] = [dx,dy]
    return distancematrix, directionmatrix
    
    
    
    
    
Eulerpropstep(system)