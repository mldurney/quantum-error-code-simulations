import glob
import os
import sys
from average_runs import rw_averages


def find_unaveraged_csv(directory):

    old_dir = os.getcwd()
    magnetization_dir = os.path.join(directory, 'magnetizations')

    if os.path.isdir(magnetization_dir):
        directory = magnetization_dir

    os.chdir(directory)
    filenames = []

    for filename in glob.glob('[!avg]*.csv'):
        filenames.append(str(os.path.join(directory, filename)))

    os.chdir(old_dir)
    return filenames

# end find_unaveraged_csv


def rw_averages_dir(directory):

    filenames = find_unaveraged_csv(directory)

    for filename in filenames:
        rw_averages(filename)

# end rw_averages_dir


def main():

    if len(sys.argv) == 2:
        directory = os.path.join(os.getcwd(), sys.argv[1])

        if not os.path.isdir(sys.argv[1]):
            print('Directory does not exist! Need magnetization directory')
            print('Usage: ' + sys.argv[0] + ' magnetizations_file_dir\n')
            sys.exit(1)

    elif len(sys.argv) == 1:
        directory = os.getcwd()

    else:
        print('Usage: ' + sys.argv[0] + ' magnetizations_file_dir\n')
        sys.exit(1)

    rw_averages_dir(directory)

# end main


if __name__ == '__main__':
    main()
