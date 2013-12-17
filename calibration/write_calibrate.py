#!/usr/bin/env python

import ConfigParser
import re
import os.path, sys

# check that the data file exists
if not os.path.exists('data/obs.csv'):
    sys.exit('data/obs.csv is missing; calibrate.m will fail')

conf = ConfigParser.ConfigParser()
conf.read('../lake.cfg')

# split the parameter names and values
param_names, param_values = zip(*conf.items('Simulation parameters'))
param_names_list = ", ".join(param_names)
param_values_list = ", ".join(param_values)
n_params = len(param_names)

# write the base calibration script
with open('calibrate.m.template', 'r') as f:
    template = f.read()

script_content = template
script_content = re.sub('__VARIABLE_NAMES_HERE__', param_names_list, script_content)
script_content = re.sub('__VARIABLE_VALUES_HERE__', param_values_list, script_content)
script_content = re.sub('__N_VARIABLES__', str(n_params), script_content)

with open('calibrate.m', 'w') as f:
    f.write(script_content)

print 'calibrate.m written using current values from ../lake.cfg'
