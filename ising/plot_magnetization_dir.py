import sys
import glob
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import AutoMinorLocator


def find_averaged_csv(directory):

    old_dir = os.getcwd()
    magnetization_dir = os.path.join(directory, 'magnetizations')

    if os.path.isdir(magnetization_dir):
        directory = magnetization_dir

    os.chdir(directory)
    filenames = []

    for filename in glob.glob('avg*.csv'):
        filenames.append(str(os.path.join(directory, filename)))

    os.chdir(old_dir)
    return filenames

# end find_averaged_csv


def make_plot_directory(directory):

    main_dir = (os.path.dirname(directory)
                if directory.endswith('magnetizations') else directory)

    plot_dir = os.path.join(main_dir, 'plots')

    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    return plot_dir

# end make_plot_directory


def plot_magnetization(filename, directory):

    old_dir = os.getcwd()
    os.chdir(directory)

    basename = os.path.basename(filename).rsplit('.', 1)[0]
    fig_filename = os.path.join(directory, basename + '.png')

    data = pd.read_csv(filename)
    temp, magn = data.columns.tolist()

    size = 5
    colors = cm.rainbow(np.linspace(0, 1, len(data[temp])))
    x_minor_locator = AutoMinorLocator(10)
    y_minor_locator = AutoMinorLocator(10)

    fig, ax = plt.subplots()
    ax.scatter(data[temp], data[magn], size, colors)
    ax.grid(color='.75', linestyle='-', linewidth=.5)
    ax.grid(which='minor', color='.01', linestyle='-', linewidth=.05)
    ax.xaxis.set_minor_locator(x_minor_locator)
    ax.yaxis.set_minor_locator(y_minor_locator)
    ax.set_title('<M> v. T: ' + basename)
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Magnetization <M>')

    fig.savefig(fig_filename, dpi=300)

    os.chdir(old_dir)



# end plot_magnetization


def main():

    if len(sys.argv) == 2:
        magnetization_directory = os.path.join(os.getcwd(), sys.argv[1])

        if not os.path.isdir(sys.argv[1]):
            print('Directory does not exist! Need magnetization directory')
            print('Usage: ' + sys.argv[0] + ' magnetizations_file_dir\n')
            sys.exit(1)

    elif len(sys.argv) == 1:
        magnetization_directory = os.getcwd()

    else:
        print('Usage: ' + sys.argv[0] + ' magnetizations_file_dir\n')
        sys.exit(1)

    filenames = find_averaged_csv(magnetization_directory)
    plot_directory = make_plot_directory(magnetization_directory)

    for filename in filenames:
        plot_magnetization(filename, plot_directory)

# end main


if __name__ == '__main__':
    main()
