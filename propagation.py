import numpy as np
from simulation import *

K = 2 #Spring constant
NEIGHBOUR_CUTOFF = 2.7 * MEAN_RADIUS
K_SELF = 1 #for now
K_BOUNDARY = 1
K_REPULSION = 1

def Eulerpropstep(system):#positions,velocities, more?):
    distances,directions = get_distances(system)    #get distances between particles and the direction of that vector
    neighbours,k,l = get_neighbours(system,distances)#get neighbours of each particle
    print(k,l)
    
    F = get_forces(system,distances, directions)
    #update position: position += timestep*velocity
    #check nearest neighbours and check if boundary (defined by angle > pi)
    #calculate average angle
    #calculate force: F = Fself + Fboundary + Frepulsion
    #update velocity: velocity += timestep*force
    #print(F[3])
    return #positions,velocities,angle?,forces
    
def get_forces(system,distances,directions):
    forcematrix = np.zeros(shape=(N_PARTICLES,N_PARTICLES,2))
    force_self, force_boundary, force_repulsion = np.zeros[2]
    for i in range(N_PARTICLES):
        for j in range(N_PARTICLES):
            angle = system[i,COLUMN_REVERSE_MAPPING['orientation']]
            orientation = [np.cos(angle),np.sin(angle)]
            #calculate self-propulsion force
            vx = system[i,COLUMN_REVERSE_MAPPING['vx']] 
            vy = system[i,COLUMN_REVERSE_MAPPING['vy']] 
            force_self = system[i,COLUMN_REVERSE_MAPPING['r']]*K_SELF*[vx,vy]
            #calculate boundary force
            if np.greater_equal(system[i,COLUMN_REVERSE_MAPPING['angle_boundary']], 180):    #if particle is part of the boundary
                force_boundary = K_BOUNDARY*(system[i,COLUMN_REVERSE_MAPPING['angle_boundary']]-180)*orientation
            #calculate repulsion force
            rsum = system[i,COLUMN_REVERSE_MAPPING['r']] + (    
                                system[j,COLUMN_REVERSE_MAPPING['r']]) #r1 + r2
            dr = distances[i,j] #distance between particle centres
            if i!=j and dr <= rsum:
                force_repulsion = -K_REPULSION* (rsum/dr - 1)*directions[i,j]
            else: force_repulsion = 0
            #calculate total force on particle
            forcematrix[i,j] = (force_self + force_boundary + force_repulsion)*system
    return forcematrix
    
def get_neighbours(system,distances):
    """Determines neighbours particles by comparing if 
        the center to center distance is <= 2.7 * MEAN_RADIUS """
        
    neighbours = np.zeros(shape=(N_PARTICLES,N_PARTICLES), dtype = bool)
    for i in range(N_PARTICLES):
        for j in range(N_PARTICLES):
            neighbours[i,j] = np.less_equal(distances[i,j],NEIGHBOUR_CUTOFF)
            if neighbours[i,j] == True:
                if np.logical_not(np.not_equal(i,j)) == False:
                    k = np.arange(i)
                    l = np.arange(j)
                
                    
    return neighbours,k,l

def Boundarycheck(system, k, l):

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
                dx = system[j,COLUMN_REVERSE_MAPPING['x']] - (
                                        system[i,COLUMN_REVERSE_MAPPING['x']])
                dy = system[j,COLUMN_REVERSE_MAPPING['y']] - (
                                        system[i,COLUMN_REVERSE_MAPPING['y']])
                distancematrix[i,j] = np.sqrt(np.power(dx,2) + np.power(dy,2))
                
                dx = dx / distancematrix[i,j]
                dy = dy / distancematrix[i,j]
                directionmatrix[i,j] = [dx,dy]
                directionmatrix[j,i] = directionmatrix[i,j] * -1

    return distancematrix, directionmatrix
    
    
    
    
    
F = get_forces(system,distances, directions)
print(F)