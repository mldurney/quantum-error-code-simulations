import os.path
from subprocess import call
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

if not os.path.isfile('average_magnetizations.csv'):
    size = '30'
    t = '.1'
    dt = '.1'
    n = '50'
    trials = '1000'
    mode = 'r'

    call('make')
    call(['./ising_simulation', size, t, dt, n, trials, mode])
    call(['make', 'clean'])

data = pd.read_csv('average_magnetizations.csv', header=None)
(m, t) = (data[0], data[1])

plt.scatter(t, m, c='blue', label='data')
plt.savefig('average_magnetizations.png', dpi=300)
plt.show()
