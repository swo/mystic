#!/usr/bin/env python

import numpy as np
import subprocess
import ConfigParser
import cPickle as pickle
import scipy.optimize
from time import gmtime, strftime
import os
import sys

log_fn = 'log'

# check that matlab path is set
if 'MATLABPATH' not in os.environ:
    die_message = "matlab path is not set! maybe source that bash script above first..."
    with open(log_fn, 'w') as f:
        f.write(die_message)
    sys.exit(die_message)

data_fn = 'tmp'
data_fn_matlab = "'{0}'".format(data_fn)
log_fn = 'log'
#matlab_args = ['matlab', '-nojvm', '-r']
matlab_args = ['matlab', '-nodesktop', '-nosplash', '-r']

# read in the base values from the config file
conf = ConfigParser.ConfigParser()
conf.read('../lake.cfg')
params = [param for name, param in conf.items('Simulation parameters')]

# read in the expected data
expected = np.genfromtxt('data/obs.csv', delimiter=',')

def matlab_f(arg_vector):
    # write out the command to send to matlab
    args = [str(arg) for arg in arg_vector] + [data_fn_matlab]
    command = matlab_args + ['get_concs({0}); exit;'.format(",".join(args))]

    # run the command
    print "calling...", " ".join(command)
    x = raw_input()
    subprocess.call(command)

    # get the data from the output file
    dat = np.genfromtxt(data_fn, delimiter=',')
    with open('dat.pkl', 'wb') as f:
        pickle.dump(dat, f)

    # get only the columns of interest
    # right now they are oxygen (0), nitrate (2), sulfate (4)
    obs = dat[:, np.array([0, 2, 4])]
    with open('obs.pkl', 'wb') as f:
        pickle.dump(obs, f)

    # get the sum of squares of differences
    sqd = np.sum(np.square(expected - obs))
    
    log_line = " ".join(command) + "\t" + str(sqd) + "\n"
    with open(log_fn, 'a') as f:
        f.write(log_line)

    print "proceed? (n stops)"
    x = raw_input()
    if x == "n":
        sys.exit('you told me to quit')

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
