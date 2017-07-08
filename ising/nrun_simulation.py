import os.path
import sys
from subprocess import call
from multiprocessing import Process
from run_simulation import run_trial


def run_trials_linear(filename, updates):

    for n in range(updates):
        print('INITIATING TRIAL #' + n)
        run_trial(filename)
        print()

# end run_trials_parallel


def run_trials_parallel(filename, updates):

    processes = [Process(target=run_trial, args=(filename,))
                 for x in range(updates)]

    for process in processes:
        process.start()

    for process in processes:
        process.join()

# end run_trials_parallel


def main():

    if len(sys.argv) != 3:
        print('Usage: ' + sys.argv[0] + ' input_file updates(int)\n')
        sys.exit(1)

    if not os.path.isfile(sys.argv[1]):
        print('Input file does not exist! Must run generate_input.py')
        print('Usage: ' + sys.argv[0] + ' input_file updates(int)\n')
        sys.exit(1)

    filename = str(sys.argv[1])
    updates = int(sys.argv[2])

    call(['make', 'clean'])
    call('make')

    run_trials_linear(filename, updates)

    call(['make', 'clean'])

# end main


if __name__ == '__main__':
    main()
