import sys
from subprocess import call
import time

if len(sys.argv == 9):
    arguments = [str(arg) for arg in sys.argv]
    (rows, cols, t, dt, n, trials, coupling, mode) = arguments[1:]

elif len(sys.argv) != 1:
    print('Usage: ' + sys.argv[0] + ' rows(int) cols(int) '
          + 'min_temperature(float) change_temperature(float) '
          + 'num_lattices(int) coupling(int) mode(char)')
    print('Alternate: ' + sys.argv[0])
    sys.exit(1)

else:
    rows = '20'
    cols = '20'
    t = '.05'
    dt = '.05'
    n = '100'
    trials = '5000'
    coupling = '1'
    mode = 'r'

call(['make', 'clean'])
call('make')

start_time = time.time()
call(['./ising_simulation', rows, cols, t, dt, n, trials, coupling, mode])
time_elapsed = time.time() - start_time

call(['make', 'clean'])

time_logs = open('times.txt', 'a')
time_logs.write('%d\n' % time_elapsed)
time_logs.close()
