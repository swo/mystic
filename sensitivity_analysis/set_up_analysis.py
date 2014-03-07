#!/usr/bin/env python

import numpy as np, re, csv, itertools, ConfigParser

# open the configuration file
conf = ConfigParser.ConfigParser()
conf.read('../lake.cfg')

# extract items from the config file
base_params = conf.items('Simulation parameters')
multipliers = np.linspace(conf.getfloat('Sensitivity analysis', 'lower_multiplier'), conf.getfloat('Sensitivity analysis', 'upper_multiplier'), conf.getint('Sensitivity analysis', 'n_multipliers'))

# prepare the command list
n_submits = conf.getint('Sensitivity analysis', 'n_submits')
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
        rates_fn = conf.get('Sensitivity analysis', 'matlab_rates_fn_mask').format(param_i, val_i)
        concs_fn = conf.get('Sensitivity analysis', 'matlab_concs_fn_mask').format(param_i, val_i)
        these_params.append(rates_fn)
        these_params.append(concs_fn)

        # stringify the filename
        params_string = ','.join([str(x) for x in these_params])

        # make the command
        command = conf.get('Sensitivity analysis', 'command_mask').format(params_string)
        commands[command_cycler.next()].append(command)

    # write out the value map
    with open(conf.get('Sensitivity analysis', 'valmap_mask').format(param_i), 'w') as f:
        w = csv.writer(f)
        w.writerows(val_list)

# write out the command files
for submit_i, commands in commands.items():
    submit_fn = conf.get('Sensitivity analysis', 'submit_fn_mask').format(submit_i)

    with open(submit_fn, 'w') as f:
        f.write("\n".join(commands))
