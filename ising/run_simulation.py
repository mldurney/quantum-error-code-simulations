import os.path
import sys
from subprocess import call

if len(sys.argv) != 2:
    print('Usage: ' + sys.argv[0] + ' input_file\n')
    sys.exit(1)

if not os.path.isfile(sys.argv[1]):
    print('Input file does not exist! Must run generate_input.py')
    print('Usage: ' + sys.argv[0] + ' input_file\n')
    sys.exit(1)

sim_name = './ising_simulation'
call('make')

tests = open(sys.argv[1], 'r')

for line in tests:
    command = [sim_name] + [setting.strip() for setting in line.split(',')]
    call(command)

call(['make', 'clean'])
