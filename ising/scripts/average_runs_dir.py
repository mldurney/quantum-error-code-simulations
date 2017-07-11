import glob
import os
import sys
from average_runs import rw_averages

FOLDERS = ['magnetizations', 'binder_cumulants']


def find_unaveraged_csv(directory, folder):

    main_dir = os.getcwd()
    folder_dir = os.path.join(directory, folder)

    if os.path.isdir(folder_dir):
        directory = folder_dir

    os.chdir(directory)
    filenames = []

    for filename in glob.glob('*.csv'):
        if 'avg' not in filename:
            filenames.append(str(os.path.join(directory, filename)))

    os.chdir(main_dir)
    return filenames

# end find_unaveraged_csv


def rw_averages_dir(directory, folder):

    filenames = find_unaveraged_csv(directory, folder)

    for filename in filenames:
        rw_averages(filename)

# end rw_averages_dir


def main():

    if len(sys.argv) == 2:
        directory = os.path.join(os.getcwd(), sys.argv[1])

        if not os.path.isdir(sys.argv[1]):
            print('Directory does not exist! Need data directory')
            print('Usage: ' + sys.argv[0] + ' data_dir\n')
            sys.exit(1)

    elif len(sys.argv) == 1:
        directory = os.getcwd()

    else:
        print('Usage: ' + sys.argv[0] + ' data_dir\n')
        sys.exit(1)

    for folder in FOLDERS:
        rw_averages_dir(directory, folder)

# end main


if __name__ == '__main__':
    main()
