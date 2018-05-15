import numpy as np
from simulation import *

K_SELF = 1.0 #for now
K_BOUNDARY = 1.0
K_REPULSION = 1.0

linear_viscosity = 1.0
angular_viscosity = 1.0

def Eulerpropstep(system):#positions,velocities, more?):
    distances,directions = get_distances(system)    #get distances between particles and the direction of that vector
    neighbours,k,l = get_neighbours(system,distances)#get neighbours of each particle
    F = get_forces(system,distances, directions)
    T = get_torque(system,distances, directions)
    #update position: position += timestep*velocity
    #update velocity: velocity += timestep*force
    #update orientation
    return #positions,velocities,angle?,forces
    
def update_velocity(system,F):
    vnew, wnew = np.zeros(shape=(2)), np.zeros(shape=(2))
    for i in range(N_PARTICLES):
        vx = system[i,COLUMN_REVERSE_MAPPING['vx']]
        vy = system[i,COLUMN_REVERSE_MAPPING['vy']]
        r  = system[i,COLUMN_REVERSE_MAPPING['r']]
        vnew[0] = vx + linear_viscosity*F[0]/r
        vnew[1] = vy + linear_viscosity*F[1]/r
        system[i,COLUMN_REVERSE_MAPPING['vx']] = vnew[0]
        system[i,COLUMN_REVERSE_MAPPING['vy']] = vnew[1]
        
    

def get_forces(system,distances,directions):
    """Calculates the net force on each particle due to its self propulsion,
        the boundary condition and the repulsion due to other particles"""
    forcematrix = np.zeros(shape=(N_PARTICLES,2))
    for i in range(N_PARTICLES):
        a = [0.0,0.0]
        force_self, force_boundary,force_repulsion,fselfandboundary = (a,a,a,a)
        #calculate self-propulsion force
        force_self = system[i,COLUMN_REVERSE_MAPPING['r']] * K_SELF
        #calculate boundary force
        outer_angle = system[i,COLUMN_REVERSE_MAPPING['angle_boundary']]
        #if particle is part of the boundary
        if np.greater_equal(outer_angle, 180.0): 
            force_boundary = K_BOUNDARY*(outer_angle - 180.0)
        #calculate repulsion force
        repulsion = np.zeros(shape=(N_PARTICLES,2))
        for j in range(N_PARTICLES):
            ri = system[i,COLUMN_REVERSE_MAPPING['r']]
            rj = system[j,COLUMN_REVERSE_MAPPING['r']]
            rsum = ri + rj
            dr = distances[i,j] #distance between particle centres
            if i!=j and dr <= rsum:
                #force on particle i by every other particle
                repulsion = -K_REPULSION*(rsum/dr - 1) * directions[i,j]
            else: repulsion = [0,0]
            force_repulsion = np.add(force_repulsion,repulsion)
        #calculate total force on particle
        angle = system[i,COLUMN_REVERSE_MAPPING['orientation']]
        orientation = [np.cos(angle),np.sin(angle)] 
        fselfandboundary[0] = (force_self + force_boundary) * orientation[0]
        fselfandboundary[1] = (force_self + force_boundary) * orientation[1]
        forcematrix[i] = np.add(fselfandboundary, force_repulsion)
    return forcematrix

    
F = get_forces(system,distances, directions)
