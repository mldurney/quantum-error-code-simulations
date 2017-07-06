import sys
from subprocess import call
import time

if len(sys.argv) == 7:
    arguments = [str(arg) for arg in sys.argv]
    (filename, t, dt, n, updates, mode) = arguments[1:]

elif len(sys.argv) != 2:
    print('Usage: ' + sys.argv[0] + ' input_file min_temperature(float)'
          + 'change_temperature(float) num_lattices(int) updates(int)'
          + 'mode(char)')
    print('Alternate: ' + sys.argv[0] + ' input_file')
    sys.exit(1)

else:
    filename = sys.argv[1]
    t = '.1'
    dt = '.1'
    n = '100'
    updates = '10000'
    coupling = '1'
    mode = 'r'


call('make')

start_time = time.time()
call(['./isingsimulation', filename, t, dt, n, updates, mode])
time_elapsed = time.time() - start_time

call(['make', 'clean'])

time_logs = open('times.txt', 'a')
time_logs.write('%d\n' % time_elapsed)
time_logs.close()
