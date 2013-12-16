#!/usr/bin/env python

'''
Reads base values from sens.cfg and writes the matlab script that calls these values. This
script must be called before using run_interactive.m
'''

import ConfigParser
import re

conf = ConfigParser.ConfigParser()
conf.read('../lake.cfg')

# construct the string of default argument values
arg_string = "\n".join(["\t{0},\t... {1}".format(value, name) for name, value in conf.items('Simulation parameters')])

# remove the last comma to prevent matlab from complaining
arg_string = re.sub(",(?=[^,]*$)", "", arg_string)

command_string = "[time_slices, concs_history, rates_history] = lake( ...\n" + arg_string + "\n);"

with open('run_interactive_defaults.m', 'w') as f:
    f.write(command_string)
