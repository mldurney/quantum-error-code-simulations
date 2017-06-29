import sys
import os.path
import itertools
import generate_hamiltonian as gh


def make_directories(folder):

    main_dir = folder
    hamiltonian_dir = os.path.join(main_dir, 'hamiltonians')
    magnetization_dir = os.path.join(main_dir, 'magnetizations')

    if not os.path.exists(main_dir):
        os.mkdir(main_dir)

    if not os.path.exists(hamiltonian_dir):
        os.mkdir(hamiltonian_dir)

    if not os.path.exists(magnetization_dir):
        os.mkdir(magnetization_dir)

    return (main_dir, hamiltonian_dir, magnetization_dir)

# end make_directories


def main():

    if len(sys.argv) == 12:
        (shape, min_rows, max_rows, min_cols, max_cols, min_temp, change_temp,
         num_temps, updates, couplings, mode) = sys.argv[1:]

        range_rows = range(min_rows, max_rows + 1)
        range_cols = range(min_cols, max_cols + 1)

    elif len(sys.argv) == 1:
        shape = str(input('Shape of lattice (s, r, t): '))

        while (shape not in [gh.SQUARE, gh.RECTANGLE, gh.TRIANGLE]):
            shape = str(input('Invalid input. Please enter s, r, or t: '))

        min_rows = int(input('Min rows: '))
        max_rows = int(input('Max rows: '))
        row_skip = int(input('Rows to skip between simulations: '))
        min_cols = int(input('Min columns: '))
        max_cols = int(input('Max columns: '))
        col_skip = int(input('Columns to skip between simulations: '))
        min_temp = float(input('Starting temperature: '))
        max_temp = float(input('Ending temperature: '))
        change_temp = float(input('Change in temperature between tests: '))
        updates = int(input('Number of updates per temperature per test '))
        coup = input('List (delimited by space) of couplings between pairs: ')
        mode = str(input('Update mode ("a" all, "r" random): '))

        range_rows = range(min_rows, max_rows + 1, row_skip)
        range_cols = range(min_cols, max_cols + 1, col_skip)
        num_temps = int((max_temp - min_temp) / change_temp) + 1
        couplings = list(map(int, coup.split()))

    else:
        shape = gh.SQUARE
        (min_rows, max_rows) = (10, 30)
        (min_cols, max_cols) = (10, 30)
        range_rows = range(min_rows, max_rows + 1)
        range_cols = range(min_cols, max_cols + 1)

        (min_temp, change_temp, num_temps) = (1.1, .1, 20)
        updates = 2000
        couplings = [1]
        mode = 'r'


    folder = (shape +  str(min_rows) + 'x' + str(min_cols) + '-'
              + str(max_rows) + 'x' + str(max_cols) + '_' + str(min_temp)
              + 'T-' + str(change_temp) + 'dTx' + str(num_temps) + '-'
              + str(updates) + 'u-' + mode)

    tests = list(itertools.product(range_rows, range_cols, couplings))


    (main_dir, hamiltonian_dir, *_) = make_directories(folder)

    tests_file = open(os.path.join(main_dir, 'tests.csv'), 'w')

    for test in tests:
        str_rows = str(test[0])
        str_cols = str(test[1])
        str_coupling = str(test[2])

        hamiltonian_name = (shape + str_rows + 'x' + str_cols + '-'
                            + str_coupling + '.csv')
        hamiltonian_filename = os.path.join(hamiltonian_dir, hamiltonian_name)

        tests_file.write('%s,%f,%f,%d,%d,%c\n' % (hamiltonian_filename,
                         min_temp, change_temp, num_temps, updates, mode))

        if shape == gh.SQUARE:
            hamiltonian = gh.generate_square(test[2], test[0])
        elif shape == gh.RECTANGLE:
            hamiltonian = gh.generate_rectangle(test[2], test[0], test[1])
        else:
            hamiltonian = gh.generate_triangle(test[2], test[0], test[1])

        gh.write_hamiltonian(hamiltonian_filename, hamiltonian, shape,
                             str_rows, str_cols)

    tests_file.close()

# end main


if __name__ == '__main__':
    main()
