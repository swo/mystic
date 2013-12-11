#!/usr/bin/env python

'''
Use the rates_...csv files produced from the matlab runs to find
the sensitivity of different positions for each rate curve to each
of the modified parameters
'''

import csv, ConfigParser
import numpy as np
import scipy.interpolate, scipy.stats

n_vals = 5

def get_single_values(x, y):
    '''give interpolated x positions of half max's, max'''

    f = scipy.interpolate.interp1d(x, y, kind='cubic')

    new_x = np.linspace(min(x), max(x), 1000)
    new_y = f(new_x)

    # find position of maximum
    max_pos = np.argmax(new_y)
    max_x = new_x[max_pos]
    max_y = np.max(new_y)

    # find the upper and lower half-max
    if np.max(new_y) - np.min(new_y) == 0:
        # escape hatch
        return [0, 0, 0]
    else:
        above_hm = new_x[new_y > (np.max(new_y) - np.min(new_y)) / 2]
        uhm = np.max(above_hm)
        lhm = np.min(above_hm)

        return [lhm, max_x, uhm]

slope = lambda x, y: scipy.stats.linregress(x, y)[0]
slopes = lambda x, ys: [slope(x, y) for y in ys]

# read in the rate names
conf = ConfigParser.ConfigParser()
conf.read('sens.cfg')
rate_names = [name for index, name in sorted(conf.items('Rate names'))]

# read in the parameter names
params = [param for param, base_val in conf.items('Simulation')]

# loop over parameters
out_rows = []
for param_i, param_name in enumerate(params):
    out_row = []

    # get values of the parameter
    fn = conf.get('Scripting', 'valmap_mask').format(param_i)
    with open(fn, 'r') as f:
        r = csv.reader(f)
        param_vals = np.array([row[1] for row in r], dtype=float)

    # loop over values of the parameter
    # initialize single values output
    single_values = np.empty((len(rate_names), len(param_vals), 3))
    for val_i, param_val in enumerate(param_vals):
        # get the data
        fn = conf.get('Scripting', 'analysis_data_fn_mask').format(param_i, val_i)
        dat = np.genfromtxt(fn, delimiter=',')
        x = range(dat.shape[0])

        # loop over each process
        for rate_i, rate_name in enumerate(rate_names):
            y = dat[:, rate_i]
            single_values[rate_i, val_i, :] = np.array(get_single_values(x, y))

    # do the linear regressions
    for rate_i, rate_name in enumerate(rate_names):
        for metric_i in [0, 1, 2]:
            #swo> here's where you pick whether to make the x values based on the actual
            # values of the parameter or just based on the percentage positions
            #my_slope = slope(param_vals, single_values[rate_i, :, metric_i])
            my_slope = slope(range(len(param_vals)), single_values[rate_i, :, metric_i])
            out_row.append(my_slope)
            print 'p={} r={} m={} slope={}'.format(param_name, rate_name, metric_i, my_slope)

    out_rows.append(out_row)

# construct the output
with open('out.csv', 'w') as f:
    w = csv.writer(f)

    header = ['rate'] + ['{0}_{1}'.format(rate_name, pos_name) for rate_name in rate_names for pos_name in ['LHM', 'max', 'UHM']]
    w.writerow(header)

    for param_name, row in zip(params, out_rows):
        new_row = [param_name] + row
        w.writerow(new_row)
