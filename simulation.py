import matplotlib.pyplot as plt
import numpy as np

from random import uniform

COLUMN_MAPPING = {
    0: 'x',
    1: 'y',
    2: 'r',
    3: 'vx',
    4: 'vy',
    5: 'v_angular',
    6: 'orientation',
    7: 'angle_boundary',
    8: 'angle_in',
    9: 'angle_delta'
}
COLUMN_REVERSE_MAPPING = {v: k for (k, v) in COLUMN_MAPPING.items()}

d = 2
N_COLUMNS = len(COLUMN_MAPPING)
N_PARTICLES = 200
LATTICE_WIDTH = 10
LATTICE_LENGTH = int(N_PARTICLES / LATTICE_WIDTH)
LATTICE_CONSTANT = 2
assert LATTICE_LENGTH * LATTICE_WIDTH == N_PARTICLES
assert LATTICE_LENGTH > LATTICE_WIDTH

MEAN_RADIUS = 1
STD_RADIUS = 1 / 10

NEIGHBOUR_CUTOFF = 2.7 * MEAN_RADIUS

K_SELF = 1.0 
K_BOUNDARY = 1.0
K_REPULSION = 100.0

TORQUE_IN = 1.0
TORQUE_NOISE = 1.0
TORQUE_ALIGN = 1.0

LINEAR_VISCOSITY = 1.0
ANGULAR_VISCOSITY = 1.0

SIMULATION_STEPS = 1000
TIME_DELTA = 1e-2
PRINT_EVERY_STEPS = 1
PLOT_EVERY_STEPS = 100

def initialize_system():
    """Initializes the system in a rectangle lattice with particles 
        slightly deviated from the exact lattice positions"""
    system = np.zeros((N_PARTICLES, N_COLUMNS), dtype = np.float64)
    length_positions = np.arange(LATTICE_LENGTH) * LATTICE_CONSTANT
    width_positions = np.arange(LATTICE_WIDTH) * LATTICE_CONSTANT 
    lattice_x, lattice_y = np.meshgrid(length_positions, width_positions)
    # initialize positions
    system[:, COLUMN_REVERSE_MAPPING['x']] = lattice_x.flatten() 
    system[:, COLUMN_REVERSE_MAPPING['y']] = lattice_y.flatten()
    
    for particle in range(N_PARTICLES):
        system[particle,COLUMN_REVERSE_MAPPING['x']] = system[
                particle,COLUMN_REVERSE_MAPPING['x']] + (np.random.rand())
        system[particle,COLUMN_REVERSE_MAPPING['y']] = system[
                particle,COLUMN_REVERSE_MAPPING['y']] + (np.random.rand())
        system[particle, COLUMN_REVERSE_MAPPING['orientation']] = (
                                                        np.random.rand() * 45)
    # initialize radius
    system[:, COLUMN_REVERSE_MAPPING['r']] = np.random.normal(
        loc=MEAN_RADIUS,
        scale=STD_RADIUS,
        size=N_PARTICLES)
    return system

def get_distances(system):
    distancematrix = np.zeros(shape=(N_PARTICLES,N_PARTICLES))
    directionmatrix = np.zeros(shape=(N_PARTICLES,N_PARTICLES,2))

    for i in range(N_PARTICLES): #for every particle
        distancematrix[i,i] = 0     #by definition #not sure its needed
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
                directionmatrix[i,j] = [dx,dy] #matrix of unitary vectors
                #Given the way you calculated them I think it's not needed
                directionmatrix[j,i] = directionmatrix[i,j] * -1 
    return distancematrix, directionmatrix


def get_neighbours(system,distances):
    """Determines neighbouring particles by comparing if 
        the center to center distance is <= 2.7 * MEAN_RADIUS """
    neighbours = np.zeros(shape=(N_PARTICLES,N_PARTICLES), dtype=bool)
    neighbours_indexes = [[] for _ in range(N_PARTICLES)]
    for i in range(N_PARTICLES):
        for j in range(N_PARTICLES):
            neighbours[i,j] = np.less_equal(distances[i,j],NEIGHBOUR_CUTOFF)
            if neighbours[i,j] == True and distances[i,j] != 0.0:
                neighbours_indexes[i].append(j)
    return neighbours_indexes

def update_angles(system, directionmatrix,neighbours_indexes):
    """Obtains the angle formed between each particle and its neighbours; 
        determines the boundary (angle_out) angle,the inner (angle_in) angle,
        and the mismatch (angle_delta) between the inner angle and the 
        orientation of the particle."""

    angles_out = [[] for _ in range(N_PARTICLES)]
    for particle in neighbours_indexes:
        for neighbour_a in range(len(particle)-1):
            neighbour_b = neighbour_a + 1
            vector_a = directionmatrix[neighbours_indexes.index(particle),
                                                   particle[neighbour_a]]
            vector_b = directionmatrix[neighbours_indexes.index(particle), 
                                                   particle[neighbour_b]]
            dot_prod = np.dot([vector_a[0],vector_a[1]],
                              [vector_b[0],vector_b[1]])
            angle_ab = np.arccos(dot_prod)
            angle_out = np.rad2deg(2*np.pi - angle_ab)
            angles_out[neighbours_indexes.index(particle)].append(angle_out)

    for particle in range(N_PARTICLES):
        system[particle, COLUMN_REVERSE_MAPPING['angle_boundary']] = (
                                    max(angles_out[particle], default=360))
        system[particle,COLUMN_REVERSE_MAPPING['angle_in']] = (
            system[particle, COLUMN_REVERSE_MAPPING['angle_boundary']] / 2.)
        system[particle, COLUMN_REVERSE_MAPPING['angle_delta']] = (
                system[particle, COLUMN_REVERSE_MAPPING['angle_in']]- (
                system[particle, COLUMN_REVERSE_MAPPING['orientation']]))

    return system, angles_out

'''
TO CHECK: parameters for F_self, F_boundary and F_repulsion
F_self ~ 1
F_boundary ~ 70
F_repulsion ~ 0.1 (or 0 if there is no overlapping)
-->
Is Fboundary too large and F_repulsion too small?

TO CHECK: parameters for torque_boundary, torque_noise, torque_align
torque_boundary ~ 60-160
torque_noise ~ 1
torque_align ~ 250-300

N.B. testingvalues.py now prints these values for each particle)
'''    
    
def get_forces(system,distances,directions):
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

def update_velocity(system,forcematrix, torque_total, time_step):
    for i in range(N_PARTICLES):
        system[i,COLUMN_REVERSE_MAPPING['vx']] += forcematrix[i][0] * time_step
        system[i,COLUMN_REVERSE_MAPPING['vy']] += forcematrix[i][1] * time_step
        system[i,COLUMN_REVERSE_MAPPING['v_angular']] += (
                                                torque_total[i] * time_step)
    return system

def update_position(system, time_step):
    for i in range(N_PARTICLES):
        vx = system[i,COLUMN_REVERSE_MAPPING['vx']]
        vy = system[i,COLUMN_REVERSE_MAPPING['vy']]
        system[i,COLUMN_REVERSE_MAPPING['x']] += vx * time_step
        system[i,COLUMN_REVERSE_MAPPING['y']] += vy * time_step
    return system

def update_orientation(system, time_step):
    for i in range(N_PARTICLES):
        w  = system[i,COLUMN_REVERSE_MAPPING['v_angular']]
        system[i,COLUMN_REVERSE_MAPPING['orientation']] += w * time_step
    return system

#---------------Observables------------------------------
    #okay so the calculation makes sense cuz cos^2+sen^2 = 1 -> 1/N_part = 0.009
    #so obviously I'm not understandig how to calculate the order parameter. 
    #According to the paper is the sum of the orientation vectors, which would 
    #give a vector not an scalar
def get_order_parameter(system):
    """Calculates the orientational order parameter. The behaviour of the 
        system can be determined from it. A high order parameter it's 
        migration while a low is jammed or rotating"""

    orientation_sum = np.abs(np.sum(
                            np.deg2rad(system[:, COLUMN_REVERSE_MAPPING['orientation']])))
    order_parameter = (1 / N_PARTICLES) * orientation_sum
    return order_parameter
#---------------- Simulation-----------------------------
def simulation_loop(system):
    """Integrates the system for a given number of steps """

    for step in range(SIMULATION_STEPS):
        time_step = TIME_DELTA
        distances, directions = get_distances(system)
        neighbours_indexes = get_neighbours(system, distances)
        angles, angles_out = update_angles(system, directions, 
                                                           neighbours_indexes)
        forces = get_forces(system, distances, directions)
        torques = get_torque(system,neighbours_indexes)
        updt_velocities = update_velocity(system, forces, torques, time_step)
        updt_position = update_position(updt_velocities, time_step)
        updt_orientation = update_orientation(updt_position, time_step)
        order_parameter = get_order_parameter(updt_orientation)
        
        if step == SIMULATION_STEPS - 1 or step % PLOT_EVERY_STEPS == 0:
            print('Step',step)

            print(order_parameter)
            plot_system(updt_orientation)
            plt.show()
    return system

# -----------Plotting--------------------------
def plot_system(system):
    """Plots the position and orientation of each particle"""
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(min(system[:,COLUMN_REVERSE_MAPPING['x']]) - MEAN_RADIUS, max(system[:,COLUMN_REVERSE_MAPPING['x']]) + MEAN_RADIUS)
    ax.set_ylim(min(system[:,COLUMN_REVERSE_MAPPING['y']]) - MEAN_RADIUS, max(system[:,COLUMN_REVERSE_MAPPING['y']]) + MEAN_RADIUS)
    for x, y, r, angle in zip(system[:, COLUMN_REVERSE_MAPPING['x']],
                       system[:, COLUMN_REVERSE_MAPPING['y']],
                       system[:, COLUMN_REVERSE_MAPPING['r']],
                       system[:,COLUMN_REVERSE_MAPPING['orientation']]):
        #print('x',x,'y', y,'r', r)
        circle = plt.Circle((x, y), radius=r, fill=False)
        ax.add_artist(circle)
    #Draw orientation arrow
        head_length = 0.05
        ax.arrow(x,y, (r-head_length) * np.cos(angle), 
                         (r-head_length) * np.sin(angle), head_width=0.05, 
                         head_length = head_length, fc='k', ec='k')

# ------------Main------------------------------
if __name__ == '__main__':
    system = initialize_system()
    plot_system(system)
    plt.show()
    simulation = simulation_loop(system)


    