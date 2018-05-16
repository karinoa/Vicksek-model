import numpy as np
from simulation import *

K_SELF = 1.0 #for now
K_BOUNDARY = 1.0
K_REPULSION = 1.0

linear_viscosity = 1.0
angular_viscosity = 1.0
time_step = 0.001


def Eulerpropstep(system):
    #print(system)
    distances,directions = get_distances(system)
    neighbours = get_neighbours(system,distances)
    F = get_forces(system,distances, directions,neighbours)
    T = get_torque(system,neighbours)
    update_velocity(system,F,T)
    update_position(system)
    update_orientation(system)
    #print(system)
    #update position: position += timestep*velocity
    #update velocity: velocity += timestep*force
    #update orientation
    #calculate average angle
    return #positions,velocities,angle?,forces
    
def update_velocity(system,F,T):
    for i in range(N_PARTICLES):
        system[i,COLUMN_REVERSE_MAPPING['vx']] += F[i][0]*time_step
        system[i,COLUMN_REVERSE_MAPPING['vy']] += F[i][1]*time_step
        system[i,COLUMN_REVERSE_MAPPING['v_angular']] += T[i]*time_step
    return

def update_position(system):
    for i in range(N_PARTICLES):
        vx = system[i,COLUMN_REVERSE_MAPPING['vx']]
        vy = system[i,COLUMN_REVERSE_MAPPING['vy']]
        system[i,COLUMN_REVERSE_MAPPING['x']] += vx*time_step
        system[i,COLUMN_REVERSE_MAPPING['y']] += vy*time_step
    return
        
def update_orientation(system):
    for i in range(N_PARTICLES):
        w  = system[i,COLUMN_REVERSE_MAPPING['v_angular']]
        system[i,COLUMN_REVERSE_MAPPING['orientation']] += w*time_step
    return
    
def get_forces(system,distances,directions,neighbours):
    """Calculates the net force on each particle due to its self propulsion,
        the boundary condition and the repulsion due to other particles"""

    forcematrix = np.zeros(shape=(N_PARTICLES,2))
    for i in range(N_PARTICLES):
        a = [0.0,0.0]
        force_self,force_boundary,force_repulsion,fselfandboundary = (a,a,a,a)
        #calculate self-propulsion force
        force_self = system[i,COLUMN_REVERSE_MAPPING['r']] * K_SELF
        #calculate boundary force
        outer_angle = system[i,COLUMN_REVERSE_MAPPING['angle_boundary']]
        #if particle is part of the boundary
        if np.greater_equal(outer_angle, 180.0): 
            force_boundary = K_BOUNDARY*(outer_angle - 180.0)
        #calculate repulsion force
        repulsion = np.zeros(shape=(N_PARTICLES,2))
        for j in neighbours[i]:
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

def get_torque(system, neighbours_indexes):
    """ Calculates the net torque on each particle due to the boundary 
        conditions, the noise(Vicksek type) and the particles trying to 
        align it's orientations"""

    torque_boundary = np.zeros(N_PARTICLES)
    torque_noise = np.zeros(N_PARTICLES)
    torque_align = np.zeros(N_PARTICLES)
    torque_total = np.zeros(N_PARTICLES)

    for particle in neighbours_indexes:
        for neighbour_a in range(len(particle)-1):
            noise = uniform(-1,1)
            torque_noise[particle] = TORQUE_NOISE * noise
            torque_align[particle] += (
                                system[neighbours_indexes.index(particle),
                                COLUMN_REVERSE_MAPPING['orientation']]-(
                                system[particle[neighbour_a],
                                COLUMN_REVERSE_MAPPING['orientation']]))

    for particle in range(N_PARTICLES):
        heavy_side = (system[particle, COLUMN_REVERSE_MAPPING[
                                                    'angle_boundary']] - 180)
        if heavy_side >= 0:
            torque_boundary[particle] = TORQUE_IN * (
                    system[particle, COLUMN_REVERSE_MAPPING['angle_delta']])
        else:
            torque_boundary[particle] = 0

        torque_total[particle] = torque_boundary[particle] + (
                                 torque_noise[particle]) + (
                                 TORQUE_ALIGN * torque_align[particle])

    return torque_total
    
Eulerpropstep(system)