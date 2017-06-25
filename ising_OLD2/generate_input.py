import sys
import itertools

if len(sys.argv) == 11:
    (min_rows, max_rows, min_cols, max_cols, min_temp, change_temp,
            num_temps, trials, couplings, mode) = sys.argv[1:]

    range_rows = range(min_rows, max_rows + 1)
    range_cols = range(min_cols, max_cols + 1)

elif len(sys.argv) == 1:
    min_rows = int(input('Min number of rows: '))
    max_rows = int(input('Max number of rows: '))
    row_skip = int(input('Rows to skip between simulations: '))
    min_cols = int(input('Min number of columns: '))
    max_cols = int(input('Max number of columns: '))
    col_skip = int(input('Columns to skip between simulations: '))
    min_temp = float(input('Starting temperature: '))
    max_temp = float(input('Ending temperature: '))
    change_temp = float(input('Change in temperature between tests: '))
    trials = int(input('Number of trials per temperature per test '))
    c = input('List (separated by single space) of couplings between pairs: ')
    mode = str(input('Update mode ("a" all indices, "r" random indices): '))

    range_rows = range(min_rows, max_rows + 1, row_skip)
    range_cols = range(min_cols, max_cols + 1, col_skip)
    num_temps = int((max_temp - min_temp) / change_temp) + 1
    couplings = list(map(int, c.split()))

else:
    (min_rows, max_rows) = (10, 30)
    (min_cols, max_cols) = (10, 30)
    range_rows = range(min_rows, max_rows + 1)
    range_cols = range(min_cols, max_cols + 1)

    (min_temp, change_temp, num_temps) = (.1, .1, 50)
    trials = 10000
    couplings = [1, -1]
    mode = 'r'

tests = list(itertools.product(range_rows, range_cols, couplings))

tests_file = open('input1.csv', 'w')

for test in tests:
    tests_file.write('%d,%d,%f,%f,%d,%d,%d,%c\n' % (test[0], test[1], min_temp,
                     change_temp, num_temps, trials, test[2], mode))

tests_file.close()
