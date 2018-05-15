import numpy as np
from simulation import *

K_SELF = 1.0 #for now
K_BOUNDARY = 1.0
K_REPULSION = 1.0

def get_forces(system,distances,directions):
    """Calculates the net force on each particle due to its self propulsion,
        the boundary condition and the repulsion due to other particles"""
    forcematrix = np.zeros(shape=(N_PARTICLES,2))
    for i in range(N_PARTICLES):
        f = [0.0,0.0]
        a,b,c,d = [],[],[],[]
        force_self, force_boundary,force_repulsion,fselfandboundary = (f,f,f,f)
        #calculate self-propulsion force
        force_self = system[i,COLUMN_REVERSE_MAPPING['r']] * K_SELF
        a.append(force_self)
        #calculate boundary force
        outer_angle = system[i,COLUMN_REVERSE_MAPPING['angle_boundary']]
        #if particle is part of the boundary
        if np.greater_equal(outer_angle, 180.0): 
            force_boundary = K_BOUNDARY*(outer_angle - 180.0)
            b.append(force_boundary)
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
        d.append(force_repulsion)
        #calculate total force on particle
        angle = system[i,COLUMN_REVERSE_MAPPING['orientation']]
        orientation = [np.cos(angle),np.sin(angle)] 
        fselfandboundary[0] = (force_self + force_boundary) * orientation[0]
        fselfandboundary[1] = (force_self + force_boundary) * orientation[1]
        c.append([fselfandboundary[0],fselfandboundary[1]])
        forcematrix[i] = np.add(fselfandboundary, force_repulsion)
    print("min, max and mean for f_self, f_boundary, both and f_repulsion are respectively:")
    print(np.min(a),np.max(a),np.average(a))
    print(np.min(b),np.max(b),np.average(b))
    print(np.min(c),np.max(c),np.average(c))
    print(np.min(d),np.max(d),np.average(d))
    return forcematrix

TORQUE_IN = 1.0
TORQUE_NOISE = 1.0
TORQUE_ALIGN = 1.0

def get_torque(system, neighbours_indexes):
    """ Calculates the net torque """
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
        if heavy_side > 0:
            torque_boundary[particle] = TORQUE_IN * (
                    system[particle, COLUMN_REVERSE_MAPPING['angle_delta']])
        else:
            torque_boundary[particle] = 0

        torque_total[particle] = torque_boundary[particle] + (
                                        torque_noise[particle]) + (
                                        TORQUE_ALIGN * torque_align[particle])
    (a,b,c) = (torque_boundary,torque_noise,TORQUE_ALIGN*torque_align) 
    print("min, max and mean for torque_boundary, torque_noise and torque_align are respectively:")
    print(np.min(a),np.max(a),np.average(a))
    print(np.min(b),np.max(b),np.average(b))
    print(np.min(c),np.max(c),np.average(c))
    return torque_total
    
F = get_forces(system,distances, directions)
print("")
T = get_torque(system,neighbours_indexes)