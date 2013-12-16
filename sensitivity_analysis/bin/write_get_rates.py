#!/usr/bin/env python

import ConfigParser

conf = ConfigParser.ConfigParser()
conf.read('../../lake.cfg')

# extract the parameter names and order
param_names = [name for name, value in conf.items('Simulation parameters')]
run_params = ", ".join(param_names)
run_sens_params = run_params + ", out_fn"

template = """function [] = get_rates({0})

[time_slices, concs_history, rates_history] = lake({1});
final_rates = rates_history(end, :, :);
dlmwrite(out_fn, final_rates);

end"""

script_content = template.format(run_sens_params, run_params)

with open('get_rates.m', 'w') as f:
    f.write(script_content)
