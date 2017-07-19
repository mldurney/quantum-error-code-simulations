import time
import os
from nrun_simulation import run_trials_parallel
import config as cf

old_dir = os.getcwd()
os.chdir(cf.MAKE_DIR)

times = [0] * 20
for _ in range(10):
    for n in range(len(times)):
        t0 = time.time()
        run_trials_parallel(
            '../s10x10-30x30_1.0T-0.1dTx21-200u-p/tests0.csv', 100, n)
        times[n] += time.time() - t0

for n in range(len(times)):
    print('Trial ' + str(n + 1) + ': ' + str(times[n] / 10))

os.chdir(old_dir)
