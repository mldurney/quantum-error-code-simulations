import os.path
import sys
import subprocess
import config as cf

FNULL = open(os.devnull, 'w')


def run_trial(filename):

    tests = open(filename, 'r')
    try:
        sim_path = cf.find_exec_path(cf.ROOT_DIR, cf.CSIM_NAME)
    except ValueError as msg:
        print(msg)
        sys.exit(1)

    for line in tests:
        args = [setting.strip() for setting in line.split(',')]
        hamiltonian_path = os.path.join(
            os.path.dirname(filename), args[0], args[1])
        del args[0:2]
        args.insert(0, hamiltonian_path)

        command = [sim_path] + args
        subprocess.call(command, stdout=FNULL, stderr=subprocess.STDOUT)

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
    #old_dir = os.getcwd()
    #os.chdir(cf.MAKE_DIR)
    #subprocess.call('make', stdout=FNULL, stderr=subprocess.STDOUT)
    run_trial(filename)
    #os.chdir(old_dir)

# end main


if __name__ == '__main__':
    main()
