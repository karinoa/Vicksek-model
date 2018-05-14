import numpy as np
from simulation import *

K = 2.0 #Spring constant
NEIGHBOUR_CUTOFF = 2.7 * MEAN_RADIUS
K_SELF = 1.0 #for now
K_BOUNDARY = 1.0
K_REPULSION = 1.0

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
    forcematrix = np.zeros(shape=(N_PARTICLES,2))
    for i in range(N_PARTICLES):
        force_self, force_boundary,force_repulsion = [0.0,0.0],[0.0,0.0],[0.0,0.0]
        #calculate self-propulsion force
        force_self = system[i,COLUMN_REVERSE_MAPPING['r']]*K_SELF
        #calculate boundary force
        angle = system[i,COLUMN_REVERSE_MAPPING['orientation']]
        orientation = [np.cos(angle),np.sin(angle)] 
        outer_angle = system[i,COLUMN_REVERSE_MAPPING['angle_boundary']]
        print(outer_angle)
        if np.greater_equal(outer_angle, 180.0):    #if particle is part of the boundary
            force_boundary = K_BOUNDARY*(outer_angle - 180.0)*orientation
        #calculate repulsion force
        repulsion = np.zeros(shape=(N_PARTICLES,2))
        for j in range(N_PARTICLES):
            ri = system[i,COLUMN_REVERSE_MAPPING['r']]
            rj = system[j,COLUMN_REVERSE_MAPPING['r']]
            rsum = ri + rj
            dr = distances[i,j] #distance between particle centres
            if i!=j and dr <= rsum:
                repulsion = -K_REPULSION*(rsum/dr - 1)*directions[i,j] #force on particle i by every other particle 
            else: repulsion = [0,0]
            force_repulsion = np.add(force_repulsion,repulsion)
        #print(force_repulsion)
            #calculate total force on particle
        #forcematrix[i] = (force_self + force_boundary) + force_repulsion
    return 0#forcematrix
    
F = get_forces(system,distances, directions)
print(F)