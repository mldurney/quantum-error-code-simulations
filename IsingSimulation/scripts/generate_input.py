'''
Generate Hamiltonians for lattice of particular shape (rectangle, square,
triangle, or square triangle), given ranges of size of lattices and integer
couplings that are the same between all neighboring pairs within a given
latytice (but may have some probability of being flipped, determiend by a
list of possible disorder values)
Generate directory for collection of lattices with subdirectories for
Hamiltonians, average magnetizations, and binder_cumulants
Output Hamiltonians as .csv file
In top directory, write tests.csv that specifies on each line a particular
Hamiltonian file that should be ran in isingsimulation along with the
temperature and update parameters for the simulation
'''

import sys
import os
import itertools
import config as cf
import generate_hamiltonian as gh


def write_hamiltonian_dir(args):
    '''
    Coordinate generation of Hamiltonian functions and associated directories
    for lattices as specified by parameters in dictionary args
    In top directory for group of lattices, write tests.csv specifying on each
    line a particular Hamiltonian file in the Hamiltonain subdirectory that
    should be rran in isingsimulation (line includes temperature and update
    parameters from running lattice)
    '''

    main_name = (args['shape'] + str(args['min_rows']) + 'x'
                 + str(args['min_cols']) + '-' + str(args['max_rows']) + 'x'
                 + str(args['max_cols']) + '_' + str(args['min_temp']) + 'T-'
                 + str(args['change_temp']) + 'dTx' + str(args['num_temps'])
                 + '-' + str(args['updates']) + 'ux' + str(args['trials'])
                 + '-' + args['mode'])
    main_dir = os.path.join(cf.ROOT_DIR, main_name)

    tests = list(itertools.product(args['range_rows'], args['range_cols'],
                                   args['couplings'], args['disorders']))

    if args['shape'] in [cf.SQUARE, cf.STRIANGLE]:
        tests = [test for test in tests if test[0] == test[1]]

    (hamiltonian_dir, *_) = make_directories(main_dir)

    tests_file = open(os.path.join(main_dir, 'tests.csv'), 'w')

    for test in tests:
        str_rows = str(test[0])
        str_cols = str(test[1])
        str_coupling = str(test[2])
        str_disorder = str(test[3])

        hamiltonian_name = (args['shape'] + str_rows + 'x' + str_cols + '_'
                            + str_coupling + '-' + str_disorder + '.csv')
        hamiltonian_filename = os.path.join(hamiltonian_dir, hamiltonian_name)

        tests_file.write('%s,%s,%f,%f,%d,%d,%d,%c\n' %
                         (os.path.basename(hamiltonian_dir), hamiltonian_name,
                          args['min_temp'], args['change_temp'],
                          args['num_temps'], args['updates'], args['trials'],
                          args['mode']))

        hamiltonian = get_hamiltonian(args['shape'], test)

        gh.write_hamiltonian(hamiltonian_filename, hamiltonian, args['shape'],
                             str_rows, str_cols)

    tests_file.close()

    return main_dir

# end write_hamiltonian_dir


def make_directories(main_dir):
    '''
    Generate directory for group of lattices in root of Ising simulation
    Create subdirectories for Hamiltonians functions and for data
    directories (i.e., magnetizations) as specified in config.py
    '''

    directories = []

    for directory in ['hamiltonians'] + cf.DIR_NAMES:
        directories.append(os.path.join(main_dir, directory))

    for directory in [main_dir] + directories:
        if not os.path.exists(directory):
            os.mkdir(directory)

    return tuple(directories)

# end make_directories


def receive_input():
    '''
    Receive input to generate input directories and files
    Return input as dictionary args
    Intended for use when generate_input.py being ran directly
    '''

    if len(sys.argv) == 14:
        args = get_cmdline_lattice()
    elif len(sys.argv) == 1:
        args = get_userinput_lattice()
    else:
        args = get_default_lattice()

    return args

# end receive_input


def get_cmdline_lattice():
    '''
    Receive input from command line arguments
    Return input as dictionary args
    '''

    args = {
        'shape': sys.argv[1],
        'min_rows': sys.argv[2],
        'max_rows': sys.argv[3],
        'min_cols': sys.argv[4],
        'max_cols': sys.argv[5],
        'min_temp': sys.argv[6],
        'change_temp': sys.argv[7],
        'num_temps': sys.argv[8],
        'updates': sys.argv[9],
        'trials': sys.argv[10],
        'couplings': sys.argv[11],
        'disorders': sys.argv[12],
        'mode': sys.argv[13]
    }

    args['range_rows'] = range(args['min_rows'], args['max_rows'] + 1)
    args['range_cols'] = range(args['min_cols'], args['max_cols'] + 1)

    return args

# end get_cmdline_lattice


def get_userinput_lattice():
    '''
    Receive input through command line interactions with user
    Return input as dictionary args
    '''

    args = {}

    args['shape'] = str(input('Shape of lattice (s, r, t, or v): '))

    while args['shape'] not in cf.SHAPES:
        args['shape'] = str(input('Invalid input. Enter s, r, t, or v: '))

    args['min_rows'] = int(input('Min rows: '))
    args['max_rows'] = int(input('Max rows: '))
    row_skip = int(input('Rows to skip between simulations: '))

    if args['shape'] in [cf.SQUARE, cf.STRIANGLE]:
        args['min_cols'] = args['min_rows']
        args['max_cols'] = args['max_rows']
        col_skip = row_skip

    else:
        args['min_cols'] = int(input('Min columns: '))
        args['max_cols'] = int(input('Max columns: '))
        col_skip = int(input('Columns to skip between simulations: '))

    args['min_temp'] = float(input('Starting temperature: '))
    max_temp = float(input('Ending temperature: '))
    args['change_temp'] = float(input('Change in temperature between tests: '))
    args['updates'] = int(
        input('Number of updates per temperature per test: '))
    args['trials'] = int(input('Number of seperate trials per lattice: '))
    couplings = input('List (delimited by space) of couplings between pairs: ')
    disorders = input('List of disorder percentages [0, 100]: ')
    args['mode'] = str(
        input('Update mode ("a" all, "p" psuedo, "r" random): '))

    args['range_rows'] = range(args['min_rows'], args['max_rows'] + 1,
                               row_skip)
    args['range_cols'] = range(args['min_cols'], args['max_cols'] + 1,
                               col_skip)
    args['num_temps'] = int((max_temp - args['min_temp'])
                            / args['change_temp']) + 1
    args['couplings'] = list(map(int, couplings.split()))
    args['disorders'] = list(map(int, disorders.split()))

    return args

# end get_userinput_lattice


def get_default_lattice():
    '''
    Generate default input
    Triggered when >= 1 command line arguments given at start of program but
    the number of arguments does not match the exact number of arguments
    needed to initialize all input values from command line
    Return input as dictionary args
    '''

    args = {
        'shape': 's',
        'min_rows': 10,
        'max_rows': 30,
        'min_temp': 1.1,
        'change_temp': .1,
        'num_temps': 20,
        'updates': 100,
        'trials': 100,
        'couplings': [1],
        'disorders': [0],
        'mode': 'r'
    }

    args['min_cols'] = args['min_rows']
    args['max_cols'] = args['max_rows']
    args['range_rows'] = range(args['min_rows'], args['max_rows'] + 1)
    args['range_cols'] = range(args['min_cols'], args['max_cols'] + 1)

    return args

# end get_default_lattice


def get_hamiltonian(shape, test):
    '''
    Call to generate appropriate hamiltonian based upon lattice shape
    and using parameters for coupling, disorder, rows, and columns as
    given in test line from tests.csv
    '''

    if shape == cf.SQUARE:
        hamiltonian = gh.generate_square(test[2], test[3], test[0])
    elif shape == cf.RECTANGLE:
        hamiltonian = gh.generate_rectangle(test[2], test[3], test[0],
                                            test[1])
    elif shape == cf.TRIANGLE:
        hamiltonian = gh.generate_triangle(test[2], test[3], test[0],
                                           test[1])
    else:
        hamiltonian = gh.generate_striangle(test[2], test[3], test[0])

    return hamiltonian

# end get_hamiltonian


if __name__ == '__main__':
    write_hamiltonian_dir(receive_input())
