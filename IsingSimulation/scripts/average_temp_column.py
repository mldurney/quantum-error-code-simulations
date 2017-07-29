import os
import sys
import re
import pandas as pd
import config as cf
from itertools import groupby
from average_runs_dir import find_unaveraged_csv

COLUMNS = ['chi0_re', 'chi0_im', 'chiq_re', 'chiq_im']


def get_lattice_name(filename):

    slash_index = filename.rfind(cf.SLASH)
    csv_index = filename.rfind('.csv')

    lattice_name = filename[slash_index + 1: csv_index]
    lattice_name = lattice_name[next(i for i, j in list(
        enumerate(lattice_name)) if j.isalpha()):]

    return lattice_name

# end get_lattice_name


def get_lattice_temperature(filename):

    slash_index = filename.rfind(cf.SLASH)
    name = filename[slash_index + 1:]
    shape_index = next(i for i, j in list(
        enumerate(name)) if j.isalpha())

    return float(name[:shape_index])

# end get_lattice_temperature


def rw_average_columns(key, data, column, directory):

    averages = {}
    for lattice in data:
        averages[lattice[0]] = lattice[1][column].mean()

    name = column + "-" + key + ".csv"
    average_filename = os.path.join(directory, name)
    df = pd.DataFrame(
        [[key, value] for key, value in averages.items()], columns=['temperature', column])
    df.to_csv(average_filename, index=False)

# end rw_average_columns


def main():

    if len(sys.argv) == 2:
        directory = os.path.join(os.getcwd(), sys.argv[1])

        if not os.path.isdir(sys.argv[1]):
            print('Directory does not exist! Need data directory')
            print('Usage: ' + sys.argv[0] + ' data_dir\n')
            sys.exit(1)

    else:
        print('Usage: ' + sys.argv[0] + ' data_dir\n')
        sys.exit(1)

    average_directory = os.path.join(directory, 'temp_averages')
    if not os.path.exists(average_directory):
        os.mkdir(average_directory)

    filenames = find_unaveraged_csv(directory, 'temp')
    filenames.sort(key=get_lattice_name)

    data = []

    for k, g in groupby(filenames, get_lattice_name):
        data = [(get_lattice_temperature(filename), pd.read_csv(filename))
                for filename in g]
        for column in COLUMNS:
            rw_average_columns(k, data, column, average_directory)

# end main


if __name__ == '__main__':
    main()
