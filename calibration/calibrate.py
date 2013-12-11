#!/usr/bin/env python

import numpy as np
import subprocess
import ConfigParser
import cPickle as pickle
import scipy.optimize
from time import gmtime, strftime

data_fn = 'tmp'
data_fn_matlab = "'{0}'".format(data_fn)
log_fn = 'log'
#matlab_args = ['matlab', '-nojvm', '-r']
matlab_args = ['matlab', '-nodesktop', '-nosplash', '-r']

# read in the base values from the config file
conf = ConfigParser.ConfigParser()
conf.read('sens.cfg')

params = [param for name, param in conf.items('Simulation')]

# read in the expected data
expected = np.genfromtxt('data/obs.csv', delimiter=',')

def matlab_f(arg_vector):
    # write out the command to send to matlab
    args = [str(arg) for arg in arg_vector] + [data_fn_matlab]
    command = matlab_args + ["run({0}); exit;".format(",".join(args))]

    # run the command
    subprocess.call(command)

    # get the data from the output file
    dat = np.genfromtxt(data_fn)

    # get the sum of squares of differences
    sqd = np.sum(np.square(expected - dat))
    
    log_line = " ".join(command) + "\t" + str(sqd)
    with open(log_fn, 'a') as f:
        f.write(log_line)

    return sqd

# clean out the log file
with open(log_fn, 'w') as f:
    time_str = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
    f.write("Calibration begun at " + time_str + "\n")

# prepare the base values as a numpy array; do the optimization
initial_values = np.array(params)

res = scipy.optimize.minimize(matlab_f, initial_values, method='Nelder-Mead')

# save the output before doing any analysis
with open('result.pkl', 'wb') as f:
    pickle.dump(res, f)

# wrap up
with open(log_fn, 'a') as f:
    time_str = strftime("%a, %d %b %Y %H:%M:%S +0000", gmtime())
    f.write("Calibration completed at " + time_str + "\n")
