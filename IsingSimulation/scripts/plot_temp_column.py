import os
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import AutoMinorLocator


def manage_plotting(directory, folder):

    filenames = find_csv(directory, folder)
    plot_directory = make_plot_directory(directory, 'temp_plots')

    for filename in filenames:
        plot_temp_data(filename, plot_directory)

# end manage_plotting


def find_csv(directory, folder):

    main_dir = os.getcwd()
    os.chdir(directory)
    folder_dir = os.path.join(os.getcwd(), folder)

    if os.path.isdir(folder_dir):
        directory = folder_dir

    os.chdir(directory)
    filenames = []

    for filename in glob.glob('*.csv'):
        filenames.append(str(os.path.join(directory, filename)))

    os.chdir(main_dir)
    return filenames

# end find_csv


def make_plot_directory(directory, folder):

    main_dir = (os.path.dirname(directory)
                if directory.endswith(folder) else directory)

    plot_dir = os.path.join(main_dir, folder)

    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    return plot_dir

# end make_plot_directory


def plot_temp_data(filename, directory):

    old_dir = os.getcwd()
    os.chdir(directory)

    plotname = os.path.basename(filename).rsplit('.', 1)[0]
    data = pd.read_csv(filename)
    temp, results = data.columns.tolist()

    size = 5
    x_minor_locator = AutoMinorLocator(10)
    y_minor_locator = AutoMinorLocator(10)
    colors = np.array([])
    c = cm.rainbow(np.linspace(0, 1, len(data[temp])))

    fig, ax = plt.subplots()
    ax.scatter(data[temp], data[results], size, c, label=plotname)

    ax.grid(color='.75', linestyle='-', linewidth=.5)
    ax.grid(which='minor', color='.01', linestyle='-', linewidth=.05)
    ax.xaxis.set_minor_locator(x_minor_locator)
    ax.yaxis.set_minor_locator(y_minor_locator)

    ax.set_title('Data v. T: ' + plotname)
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Data')
    plt.legend(loc=0)

    fig_filename = plotname + '.png'
    fig.savefig(fig_filename, dpi=300)

    os.chdir(old_dir)

# end plot_temp_data


def main():

    if len(sys.argv) == 2:
        directory = os.path.join(os.getcwd(), sys.argv[1])

    else:
        print('Usage: ' + sys.argv[0] + ' data_dir mode')
        sys.exit(1)

    if not os.path.isdir(sys.argv[1]):
        print('Directory does not exist! Need data directory')
        print('Usage: ' + sys.argv[0] + ' data_dir')
        sys.exit(1)

    manage_plotting(sys.argv[1], 'temp_averages')

# end main()


if __name__ == '__main__':
    main()
