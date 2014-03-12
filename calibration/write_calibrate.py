#!/usr/bin/env python

import ConfigParser
import re
import os, os.path, sys

class Param(object):
    def __init__(self, name, base_value, varied=False):
        self.name = name
        self.base_value = base_value
        self.varied = varied

    def __str__(self):
        if self.varied:
            return self.name
        else:
            return self.base_value

# check that the data file exists
if not os.path.exists('data/obs.csv'):
    sys.exit('data/obs.csv is missing; calibrate.m will fail')

conf = ConfigParser.ConfigParser()
conf.read('../lake.cfg')

# check that all the parameters listed as varied are actual parameters in the simulation
calibration_names = [name for name, value in conf.items('Calibration parameters')]
simulation_names = [name for name, value in conf.items('Simulation parameters')]
missing_names = set(calibration_names) - set(simulation_names)

if missing_names:
    raise RuntimeError('some calibration parameters are not simulation parameters: {0}'.format(missing_names))

# get the information about the varied values
varied_bounds = {name: [float(x) for x in value.split(',')] for name, value in conf.items('Calibration parameters')}

# split the parameter names and values. divide them up in to varied and fixed
params = [Param(name, value, varied=(name in varied_bounds)) for name, value in conf.items('Simulation parameters')]
fixed = [param for param in params if not param.varied]
varied = [param for param in params if param.varied]

# put the bounds into the varied containers
for param in varied:
    param.lb, param.ub = varied_bounds[param.name]

# prepare all the stuff to be shoved into the template
initial_varied_list = ", ".join([str(param.base_value) for param in varied])
varied_names = ", ".join([param.name for param in varied])
n_varied = str(len(varied))
call_list = ", ".join([str(param) for param in params])
lbs = ", ".join([str(param.lb) for param in varied])
ubs = ", ".join([str(param.ub) for param in varied])
method = "'{0}'".format(conf.get('Calibration settings', 'Method'))

# write the base calibration script
with open('calibrate.m.template', 'r') as f:
    template = f.read()

script_content = template
for pattern, replacement in [('__METHOD__', method), ('__INITIAL_VARIED_VALUES__', initial_varied_list), ('__VARIED_NAMES__', varied_names), ('__N_VARIED__', n_varied), ('__CALL_LIST__', call_list), ('__LBS__', lbs), ('__UBS__', ubs)]:
    script_content = re.sub(pattern, replacement, script_content)

with open('calibrate.m', 'w') as f:
    f.write(script_content)

print 'calibrate.m written using current values from ../lake.cfg'
