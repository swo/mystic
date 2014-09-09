# Mystic

An inferred biomass biogeochemical model of Mystic Lake. 

- by: Scott W Olesen (swo@mit.edu)
- development page: http://github.com/swo/mystic
- main page: http://almlab.mit.edu/mystic.html

## Dependencies

- Matlab (developed with version 8)
- Python (developed with version 2.7)

## Installation
    $ git clone http://github.com/swo/mystic.git

## File structure
- analysis: Tools for analyzing divergences between rates at each depth, diversity of rates, and the profiles for rates and concs.
- bin: Main matlab scripts
- calibration: Tools for calibrating the model parameters using the measured rates and/or OTU data.
- interactive: A script for running the model with easily-adjustable parameters and displaying the output.
- sensitivity\_analysis: Tools for analyzing the sensitivity of the rate profiles to the values for each of the input parameters.
- lake.cfg: Configuration for different tools, with lists of the names for processes and some default values from which to generate scripts for the various tools.
- set\_matlab\_path.sh: Convenience tool for setting the matlab path when running scripts from the terminal.

## Usage

### Interactive
To run the model interactively,
- Update `lake.cfg` with the desired parameters.
- Under interactive, use `write_default_values_script.py` to create `run_interactive_defaults.m`, which supplies the values in `lake.cfg` to `interactive.m`.
- In Matlab, run `run_interactive.m`. This requires the `bin/` folder to be on Matlab's path. If running in the terminal, you can use the `set_matlab_path.py` script to do this.
- The output of the simulation is stored in `concs_history` and `rates_history` variables. You can write them to files using `bin/write_data_to_files.m`.

### Analysis
To make a plot of the Jensen-Shannon divergences between the composition of rates at each depth, run the `analysis/divergenes/pipeline.sh`.

To compare the simulated Shannon diversity of composition of rates with some observed diversity of rates, run `analysis/diversity/get_diversity_depth.m`.

To make timecourses and end-time plots for all the rates and chemical concentrations, run `analysis/profiles/pipeline.sh`.


