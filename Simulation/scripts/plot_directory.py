import sys
import glob
import os
import re
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import config as cf
from matplotlib.ticker import AutoMinorLocator

plt.switch_backend('agg')

MODES = [SINGLE, COUPLING, DISORDER, SIZE, ALL] = ['s', 'c', 'd', 'z', 'a']
YMAX = 0
YMIN = 0
DEFAULT = -1


def manage_plotting(directory, mode, clean, min_temperature, max_temperature):

    for folder in cf.DIR_NAMES:
        filenames = find_averaged_csv(directory, folder)
        plot_directory = make_plot_directory(directory, folder, mode, clean)
        groups = find_filename_groups(filenames[:], mode)

        for group in groups:
            plot_results(group, plot_directory, folder, mode,
                         clean, min_temperature, max_temperature)

 # end manage_plotting


def find_averaged_csv(directory, folder):

    main_dir = os.getcwd()
    folder_dir = os.path.join(directory, folder)

    if os.path.isdir(folder_dir):
        directory = folder_dir

    os.chdir(directory)
    filenames = []

    for filename in glob.glob('*.csv'):
        filenames.append(str(os.path.join(directory, filename)))

    os.chdir(main_dir)
    return filenames

# end find_averaged_csv


def make_plot_directory(directory, folder, mode, clean):

    if mode == SINGLE:
        ending = ''
    else:
        ending = '-' + mode

    if clean:
        ending += '_clean'

    main_dir = (os.path.dirname(directory)
                if directory.endswith(folder) else directory)

    plot_dir = os.path.join(main_dir, folder + '_plots' + ending)

    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    return plot_dir

# end make_plot_directory


def find_filename_groups(filenames, mode):

    groups = []

    while filenames:
        group = []
        match = filenames[0]

        if mode == COUPLING:
            for filename in filenames:
                if (re.sub('(?<=_)\d+', 'x', match) ==
                        re.sub('(?<=_)\d+', 'x', filename)):
                    group.append(filename)

        elif mode == DISORDER:
            for filename in filenames:
                if (re.sub('(?<=-)\d+', 'x', match) ==
                        re.sub('(?<=-)\d+', 'x', filename)):
                    group.append(filename)

        elif mode == SIZE:
            for filename in filenames:
                if (re.sub('\d+x\d+', 'x', match) ==
                        re.sub('\d+x\d+', 'x', filename)):
                    group.append(filename)

        else:
            group.append(match)

        for member in group:
            filenames.remove(member)

        groups.append(group)

    return groups

# end find_filenane_groups


def plot_results(filenames, directory, folder, mode, clean, min_temperature, max_temperature):

    old_dir = os.getcwd()
    os.chdir(directory)

    size = 5
    x_minor_locator = AutoMinorLocator(10)
    y_minor_locator = AutoMinorLocator(10)
    colors = np.array([])

    fig, ax = plt.subplots()

    for i, filename in enumerate(filenames):
        basename = os.path.basename(filename).rsplit('.', 1)[0]

        print('Plotting ' + folder + cf.SLASH + basename)

        data = pd.read_csv(filename)
        temp, results = data.columns.tolist()

        if clean:
            mask = data[results] < data[results][0] * 1.1
            data = data[mask]

        if min_temperature is not DEFAULT and max_temperature is not DEFAULT:
            data = data[data[temp] >= min_temperature]
            data = data[data[temp] <= max_temperature]

        if mode == SINGLE:
            c = cm.rainbow(np.linspace(0, 1, len(data[temp])))
        else:
            if colors.size == 0:
                colors = cm.rainbow(np.linspace(0, 1, len(filenames)))
            c = colors[i]

        ax.scatter(data[temp], data[results], size, c, label=basename)

    ax.grid(color='.75', linestyle='-', linewidth=.5)
    ax.grid(which='minor', color='.01', linestyle='-', linewidth=.05)

    results_min = min(data[results])
    results_max = max(data[results])

    YMIN = results_min - results_max * .1
    YMAX = results_max * 1.1

    curr_min, curr_max = ax.get_ylim()
    if YMIN < curr_min:
        YMIN = curr_min
    if YMAX > curr_max:
        YMAX = curr_max

    ax.set_ylim([YMIN, YMAX])

    ax.xaxis.set_minor_locator(x_minor_locator)
    ax.yaxis.set_minor_locator(y_minor_locator)

    if mode == SINGLE:
        plotname = basename
    elif mode == COUPLING:
        plotname = re.sub('(?<=_)\d+', 'x', basename)
    elif mode == DISORDER:
        plotname = re.sub('(?<=-)\d+', 'x', basename)
    elif mode == SIZE:
        plotname = re.sub('\d+x\d+', 'x', basename)

    if folder == 'magnetizations':
        title = '<M> v. T: ' + plotname
        ylabel = 'Magnetization <M>'
    elif folder == 'binder_cumulants':
        title = 'U**4 v. T: ' + plotname
        ylabel = 'Binder Cumulant (U**4)'
    elif folder == 'correlation_functions':
        title = 'E_m v. T: ' + plotname
        ylabel = 'Correlation Function (E_m / L)'
    else:
        title = 'Data v. T: ' + plotname
        ylabel = 'Data'

    ax.set_title(title)
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel(ylabel)
    plt.legend(loc=0)

    fig_filename = os.path.join(directory, plotname + '.png')
    fig.savefig(fig_filename, dpi=300)

    os.chdir(old_dir)

# end plot_results


def main():

    parser = argparse.ArgumentParser(
        description='Plot results against temperature for Ising simulation data')

    parser.add_argument('directory', metavar='dir', nargs='?', default=os.getcwd(
    ), help='Main directory for specific simulation run')
    parser.add_argument('-s', '--single', action='store_true',
                        help='Plot all data separately')
    parser.add_argument('-c', '--coupling', action='store_true',
                        help='Plot data by like coupling')
    parser.add_argument('-d', '--disorder', action='store_true',
                        help='Plot data by like disorder')
    parser.add_argument('-z', '--size', action='store_true',
                        help='Plot data by like size')
    parser.add_argument('-a', '--all', action='store_true',
                        help='Plot data by all possible modes')
    parser.add_argument('-x', '--clean', action='store_true',
                        help='Remove extreme data points to reduce plot range')
    parser.add_argument('-t', '--temperatures', metavar='T', nargs=2, type=float, default=[DEFAULT, DEFAULT],
                        help='Set min and max [min max] temperature ranges')

    args = parser.parse_args()
    mode = None

    directory = os.path.join(os.getcwd(), args.directory)
    clean = True if args.clean else False
    (min_temperature, max_temperature) = args.temperatures

    if args.single:
        mode = SINGLE
        manage_plotting(directory, mode, clean,
                        min_temperature, max_temperature)
    if args.coupling:
        mode = COUPLING
        manage_plotting(directory, mode, clean,
                        min_temperature, max_temperature)
    if args.disorder:
        mode = DISORDER
        manage_plotting(directory, mode, clean,
                        min_temperature, max_temperature)
    if args.size:
        mode = SIZE
        manage_plotting(directory, mode, clean,
                        min_temperature, max_temperature)
    if args.all:
        for option in [SINGLE, COUPLING, DISORDER, SIZE]:
            mode = option
            manage_plotting(directory, mode, clean,
                            min_temperature, max_temperature)
    if mode is None:
        manage_plotting(directory, SINGLE, clean,
                        min_temperature, max_temperature)

# end main


if __name__ == '__main__':
    main()
