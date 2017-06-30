import os.path
import sys
from subprocess import call
from multiprocessing import Process
from run_simulation import run_trial


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

    processes = [Process(target=run_trial, args=(filename,))
                 for x in range(updates)]

    call('make')

    for process in processes:
        process.start()

    for process in processes:
        process.join()

    call(['make', 'clean'])

# end main


if __name__ == '__main__':
    main()
