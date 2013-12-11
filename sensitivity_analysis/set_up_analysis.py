#!/usr/bin/env python

import numpy as np, re, csv, itertools, ConfigParser

# open the configuration file
conf = ConfigParser.ConfigParser()
conf.read('sens.cfg')

# extract items from the config file
base_params = conf.items('Simulation')
multipliers = np.linspace(conf.getfloat('Analysis', 'lower_multiplier'), conf.getfloat('Analysis', 'upper_multiplier'), conf.getint('Analysis', 'n_multipliers'))

# prepare the command list
n_submits = conf.getint('Scripting', 'n_submits')
command_cycler = itertools.cycle(range(n_submits))
commands = {i: [] for i in range(n_submits)}

# modulate each of the parameters in turn
for param_i, (param_name, base_val) in enumerate(base_params):
    # when modulating this parameter, multiply by some multiplier
    val_list = []
    for val_i, multiplier in enumerate(multipliers):
        these_params = [float(base_val) for param_name, base_val in base_params]
        this_val = these_params[param_i] * multiplier
        these_params[param_i] = this_val

        val_list.append([val_i, this_val])

        # add the filename to the list of parameters
        out_fn = conf.get('Scripting', 'matlab_data_fn_mask').format(param_i, val_i)
        these_params.append(out_fn)

        # stringify the filename
        params_string = ','.join([str(x) for x in these_params])

        # make the command
        command = conf.get('Scripting', 'command_mask').format(params_string)
        commands[command_cycler.next()].append(command)

    # write out the value map
    with open(conf.get('Scripting', 'valmap_mask').format(param_i), 'w') as f:
        w = csv.writer(f)
        w.writerows(val_list)

# write out the command files
for submit_i, commands in commands.items():
    submit_fn = conf.get('Scripting', 'submit_fn_mask').format(submit_i)

    with open(submit_fn, 'w') as f:
        f.write("\n".join(commands))
