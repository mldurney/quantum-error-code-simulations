import os.path
import sys
from subprocess import call
import pandas as pd

if not os.path.isfile('input.csv'):
    sys.stderr.write('Input file does not exist! Must run generate_input.py')
    sys.exit(1)

sim_name = './ising_simulation'
call('make')

tests = open('input.csv', 'r')

for line in tests:
    command = [sim_name] + [setting.strip() for setting in line.split(',')]
    call(command)

call(['make', 'clean'])
