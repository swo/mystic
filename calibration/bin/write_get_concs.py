#!/usr/bin/env python

import ConfigParser

conf = ConfigParser.ConfigParser()
conf.read('../../lake.cfg')

# extract the parameter names and order
param_names = [name for name, value in conf.items('Simulation parameters')]
run_params = ", ".join(param_names)
run_sens_params = run_params + ", out_fn"

template = """function [] = get_concs({0})

[time_slices, concs_history, rates_history] = run({1});
final_concs = squeeze(concs_history(end, :, :));
dlmwrite(out_fn, final_concs);

end"""

script_content = template.format(run_sens_params, run_params)

with open('get_concs.m', 'w') as f:
    f.write(script_content)
