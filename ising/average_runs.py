import os.path
import sys
import pandas as pd
import numpy as np


def average_no_outliers(data, percentile1, percentile2):

    cols = ['temperature', 'magnetization']
    averages = pd.DataFrame(columns=cols)

    for temp in data.columns:
        quartile1 = np.percentile(data[temp], percentile1)
        quartile3 = np.percentile(data[temp], percentile2)
        restricted = [x for x in data[temp] if quartile1 <= x <= quartile3]
        mean = np.mean(restricted)
        std = np.std(restricted)
        no_outliers = [x for x in data[temp] if abs(x - mean) <= std]
        average = sum(no_outliers) / len(no_outliers)
        averages.loc[len(averages)] = [float(temp), average]

    return averages

# end average_no_outliers


def main():

    if len(sys.argv) != 2:
        print('Usage: ' + sys.argv[0] + ' input_file\n')
        sys.exit(1)

    if not os.path.isfile(sys.argv[1]):
        print('Input file does not exist! Must run generate_input.py')
        print('Usage: ' + sys.argv[0] + ' input_file\n')
        sys.exit(1)

    in_filename = sys.argv[1]

    data = pd.read_csv(in_filename)
    averages = average_no_outliers(data, 25, 75)

    slash_index = in_filename.rfind('\\')
    out_filename = ('avg_' + in_filename if slash_index == -1 else
                    in_filename[:slash_index + 1] + 'avg_' +
                    in_filename[slash_index + 1:])

    averages.to_csv(out_filename)

# end main


if __name__ == '__main__':
    main()
