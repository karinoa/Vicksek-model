import matplotlib.pyplot as plt
import numpy as np

COLUMN_MAPPING = {
    0: 'x',
    1: 'y',
    3: 'r',
    4: 'fx',
    5: 'fy'
}
COLUMN_REVERSE_MAPPING = {v: k for (k, v) in COLUMN_MAPPING.items()}

d = 2
N_COLUMNS = len(COLUMN_MAPPING)
N_PARTICLES = 110
LATTICE_WIDTH = 10
LATTICE_LENGTH = int(N_PARTICLES / LATTICE_WIDTH)
LATTICE_CONSTANT = 2
assert LATTICE_LENGTH * LATTICE_WIDTH == N_PARTICLES
assert LATTICE_LENGTH > LATTICE_WIDTH

MEAN_RADUIS = 1
STD_RADIUS = 1 / 10
seed = np.random

def initialize_system():
    """Initializes the system in a rectangle lattice with particles 
        slightly deviated from the exact lattice positions"""
    system = np.zeros((N_PARTICLES, N_COLUMNS))
    length_positions = np.arange(LATTICE_LENGTH) * LATTICE_CONSTANT + (np.random.rand() / STD_RADIUS)
    width_positions = np.arange(LATTICE_WIDTH) * LATTICE_CONSTANT + (np.random.rand() / STD_RADIUS)
    lattice_x, lattice_y = np.meshgrid(length_positions, width_positions)
    # initialize positions
    system[:, COLUMN_REVERSE_MAPPING['x']] = lattice_x.flatten()
    system[:, COLUMN_REVERSE_MAPPING['y']] = lattice_y.flatten()
    # initialize radius
    system[:, COLUMN_REVERSE_MAPPING['r']] = np.random.normal(
        loc=MEAN_RADUIS,
        scale=STD_RADIUS,
        size=N_PARTICLES)
    return system

def plot_system(system):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(-MEAN_RADUIS, (LATTICE_LENGTH * LATTICE_CONSTANT) +
                MEAN_RADUIS)
    ax.set_ylim(-MEAN_RADUIS, (LATTICE_WIDTH * LATTICE_CONSTANT) +
                MEAN_RADUIS)
    for x, y, r in zip(system[:, COLUMN_REVERSE_MAPPING['x']],
                       system[:, COLUMN_REVERSE_MAPPING['y']],
                       system[:, COLUMN_REVERSE_MAPPING['r']]):
        #print('x',x,'y', y,'r', r)
        circle = plt.Circle((x, y), radius=r)
        ax.add_artist(circle)

if __name__ == '__main__':
    system = initialize_system()
    plot_system(system)
    plt.show()