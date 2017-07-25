import os.path
import sys
import pandas as pd
import numpy as np
import config as cf

MIN_PERCENTILE = 25
MAX_PERCENTILE = 75


def average_no_outliers(data, percentile1, percentile2):

    cols = ['temperature', 'results']
    averages = pd.DataFrame(columns=cols)

    for temp in data.columns:
        data_no_nan = data[temp].dropna()
        quartile1 = np.percentile(data_no_nan, percentile1)
        quartile3 = np.percentile(data_no_nan, percentile2)
        restricted = [x for x in data_no_nan if quartile1 <= x <= quartile3]
        mean = np.mean(restricted)
        std = np.std(restricted)
        no_outliers = [x for x in data_no_nan if abs(x - mean) <= std]
        average = sum(no_outliers) / len(no_outliers)
        averages.loc[len(averages)] = [float(temp), average]

    return averages

# end average_no_outliers


def rw_averages(in_filename):

    # print('Entering ' + in_filename)
    data = pd.read_csv(in_filename)
    averages = average_no_outliers(data, MIN_PERCENTILE, MAX_PERCENTILE)

    slash_index = in_filename.rfind(cf.SLASH)
    out_filename = ('avg_' + in_filename if slash_index == -1 else
                    in_filename[:slash_index + 1] + 'avg_' +
                    in_filename[slash_index + 1:])

    averages.to_csv(out_filename, index=False)

# end rw_averages


def main():

    if len(sys.argv) != 2:
        print('Usage: ' + sys.argv[0] + ' input_file\n')
        sys.exit(1)

    if not os.path.isfile(sys.argv[1]):
        print('Input file does not exist! Must run generate_input.py')
        print('Usage: ' + sys.argv[0] + ' input_file\n')
        sys.exit(1)

    rw_averages(sys.argv[1])

# end main


if __name__ == '__main__':
    main()
