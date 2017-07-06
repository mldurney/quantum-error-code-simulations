import sys
import glob
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import AutoMinorLocator

plt.switch_backend('agg')

FOLDERS = ['magnetizations', 'binder_cumulants']


def find_averaged_csv(directory, folder):

    main_dir = os.getcwd()
    folder_dir = os.path.join(directory, folder)

    if os.path.isdir(folder_dir):
        directory = folder_dir

    os.chdir(directory)
    filenames = []

    for filename in glob.glob('avg*.csv'):
        filenames.append(str(os.path.join(directory, filename)))

    os.chdir(main_dir)
    return filenames

# end find_averaged_csv


def make_plot_directory(directory, folder):

    main_dir = (os.path.dirname(directory)
                if directory.endswith(folder) else directory)

    plot_dir = os.path.join(main_dir, folder + '_plots')

    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    return plot_dir

# end make_plot_directory


def plot_results(filename, directory, folder):

    old_dir = os.getcwd()
    os.chdir(directory)

    basename = os.path.basename(filename).rsplit('.', 1)[0]
    fig_filename = os.path.join(directory, basename + '.png')

    data = pd.read_csv(filename)
    temp, results = data.columns.tolist()

    size = 5
    colors = cm.rainbow(np.linspace(0, 1, len(data[temp])))
    x_minor_locator = AutoMinorLocator(10)
    y_minor_locator = AutoMinorLocator(10)

    fig, ax = plt.subplots()
    ax.scatter(data[temp], data[results], size, colors)
    ax.grid(color='.75', linestyle='-', linewidth=.5)
    ax.grid(which='minor', color='.01', linestyle='-', linewidth=.05)
    ax.xaxis.set_minor_locator(x_minor_locator)
    ax.yaxis.set_minor_locator(y_minor_locator)

    if folder == 'magnetizations':
        title = '<M> v. T: ' + basename
        ylabel = 'Magnetization <M>'
    elif folder == 'binder_cumulants':
        title = 'U**4 v. T: ' + basename
        ylabel = 'Binder Cumulant (U**4)'
    else:
        title = 'Data v. T: ' + basename
        ylabel = 'Data'

    ax.set_title(title)
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel(ylabel)

    fig.savefig(fig_filename, dpi=300)

    os.chdir(old_dir)

# end plot_magnetization


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

    folders = ['magnetizations', 'binder_cumulants']

    for folder in folders:
        filenames = find_averaged_csv(directory, folder)
        plot_directory = make_plot_directory(directory, folder)

        for filename in filenames:
            plot_results(filename, plot_directory, folder)

# end main


if __name__ == '__main__':
    main()
