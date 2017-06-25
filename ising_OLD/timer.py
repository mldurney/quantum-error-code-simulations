from subprocess import call
import time

size = '20'
t = '.1'
dt = '.1'
n = '50'
trials = '1000'
mode = 'r'

call(['make', 'clean'])
call('make')

start_time = time.time()
call(['./ising_simulation', size, t, dt, n, trials, mode])
time_elapsed = time.time() - start_time

call(['make', 'clean'])

time_logs = open('times.txt', 'a')
time_logs.write('%d\n' % time_elapsed)
time_logs.close()
