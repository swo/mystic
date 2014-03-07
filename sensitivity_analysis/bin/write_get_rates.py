#!/usr/bin/env python

import ConfigParser

conf = ConfigParser.ConfigParser()
conf.read('../../lake.cfg')

# extract the parameter names and order
param_names = [name for name, value in conf.items('Simulation parameters')]
run_params = ", ".join(param_names)
run_sens_params = run_params + ", rates_out_fn, concs_out_fn"

template = """function [] = get_rates({0})

[time_slices, concs_history, rates_history] = lake({1});

final_rates = squeeze(rates_history(end, :, :));
dlmwrite(rates_out_fn, final_rates);

final_concs = squeeze(concs_history(end, :, :));
dlmwrite(concs_out_fn, final_concs);

end"""

script_content = template.format(run_sens_params, run_params)

with open('get_rates.m', 'w') as f:
    f.write(script_content)
