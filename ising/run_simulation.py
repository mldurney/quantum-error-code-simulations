import os.path
import sys
from subprocess import call


def run_trial(filename):

    sim_name = './isingsimulation'
    tests = open(filename, 'r')

    for line in tests:
        command = [sim_name] + [setting.strip()
                                for setting in line.split(',')]
        call(command)

# end run_trial


def main():

    if len(sys.argv) != 2:
        print('Usage: ' + sys.argv[0] + ' input_file\n')
        sys.exit(1)

    if not os.path.isfile(sys.argv[1]):
        print('Input file does not exist! Must run generate_input.py')
        print('Usage: ' + sys.argv[0] + ' input_file\n')
        sys.exit(1)

    filename = str(sys.argv[1])

    call('make')

    run_trial(filename)

    call(['make', 'clean'])

# end main


if __name__ == '__main__':
    main()
