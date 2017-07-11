import sys
import os
import itertools
import config as cf
import generate_hamiltonian as gh


def write_hamiltonian_dir(args):

    main_name = (args['shape'] +  str(args['min_rows']) + 'x'
                 + str(args['min_cols']) + '-' + str(args['max_rows']) + 'x'
                 + str(args['max_cols']) + '_' + str(args['min_temp']) + 'T-'
                 + str(args['change_temp']) + 'dTx' + str(args['num_temps'])
                 + '-' + str(args['updates']) + 'u-' + args['mode'])
    main_dir = os.path.join(cf.ROOT_DIR, main_name)

    tests = list(itertools.product(args['range_rows'], args['range_cols'],
                                   args['couplings'], args['disorders']))

    if args['shape'] in [cf.SQUARE, cf.STRIANGLE]:
        tests = [test for test in tests if test[0] == test[1]]

    (hamiltonian_dir, *_) = make_directories(main_dir)
    hamiltonian_name = os.path.basename(hamiltonian_dir)

    tests_file = open(os.path.join(main_dir, 'tests.csv'), 'w')

    for test in tests:
        str_rows = str(test[0])
        str_cols = str(test[1])
        str_coupling = str(test[2])
        str_disorder = str(test[3])

        hamiltonian_fileid = (args['shape'] + str_rows + 'x' + str_cols + '_'
                              + str_coupling + '-' + str_disorder + '.csv')
        hamiltonian_filename = os.path.join(hamiltonian_dir, 
                                            hamiltonian_fileid)

        tests_file.write('%s,%s,%f,%f,%d,%d,%c\n' % 
                         (hamiltonian_name, hamiltonian_fileid, 
                          args['min_temp'], args['change_temp'], 
                          args['num_temps'], args['updates'], args['mode']))

        if args['shape'] == cf.SQUARE:
            hamiltonian = gh.generate_square(test[2], test[3], test[0])
        elif args['shape'] == cf.RECTANGLE:
            hamiltonian = gh.generate_rectangle(test[2], test[3], test[0],
                                                test[1])
        elif args['shape'] == cf.TRIANGLE:
            hamiltonian = gh.generate_triangle(test[2], test[3], test[0],
                                               test[1])
        else:
            hamiltonian = gh.generate_striangle(test[2], test[3], test[0])

        gh.write_hamiltonian(hamiltonian_filename, hamiltonian, args['shape'],
                             str_rows, str_cols)

    tests_file.close()

    return main_dir

# end write_hamiltonian_dir


def make_directories(main_dir):

    directories = []
    
    for directory in ['hamiltonians'] + cf.DIR_NAMES:
        directories.append(os.path.join(main_dir, directory))

    for directory in [main_dir] + directories:
        if not os.path.exists(directory):
            os.mkdir(directory)

    return tuple(directories)

# end make_directories


def receive_input():

    if len(sys.argv) == 13:
        args = get_cmdline_lattice()
    elif len(sys.argv) == 1:
        args = get_userinput_lattice()
    else:
        args = get_default_lattice()

    return args

# end receive_input


def get_cmdline_lattice():

    args = {
        'shape'         : sys.argv[1],
        'min_rows'      : sys.argv[2],
        'max_rows'      : sys.argv[3],
        'min_cols'      : sys.argv[4],
        'max_cols'      : sys.argv[5],
        'min_temp'      : sys.argv[6],
        'change_temp'   : sys.argv[7],
        'num_temps'     : sys.argv[8],
        'updates'       : sys.argv[9],
        'couplings'     : sys.argv[10],
        'disorders'     : sys.argv[11],
        'mode'          : sys.argv[12]
    }

    args['range_rows'] = range(args['min_rows'], args['max_rows'] + 1)
    args['range_cols'] = range(args['min_cols'], args['max_cols'] + 1)

    return args

# end get_cmdline_lattice


def get_userinput_lattice():

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
    args['updates'] = int(input('Number of updates per temperature per test: '))
    couplings = input('List (delimited by space) of couplings between pairs: ')
    disorders = input('List of disorder percentages [0, 100]: ')
    args['mode'] = str(input('Update mode ("a" all, "p" psuedo, "r" random): '))

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

    args = {
        'shape'         : 's',
        'min_rows'      : 10,
        'max_rows'      : 30,
        'min_temp'      : 1.1,
        'change_temp'   : .1,
        'num_temps'     : 20,
        'updates'       : 2000,
        'couplings'     : [1],
        'disorders'     : [0],
        'mode'          : 'r'
    }

    args['min_cols'] = args['min_rows']
    args['max_cols'] = args['max_rows']
    args['range_rows'] = range(args['min_rows'], args['max_rows'] + 1)
    args['range_cols'] = range(args['min_cols'], args['max_cols'] + 1)

    return args

# end get_default_lattice


if __name__ == '__main__':
    write_hamiltonian_dir(receive_input())
