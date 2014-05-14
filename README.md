# Structure
- analysis: Tools for analyzing divergences between rates at each depth, diversity of rates,
  and the profiles for rates and concs.
- bin: Main matlab scripts
- calibration: Tools for calibrating the model parameters using the measured rates and/or
  OTU data.
- doc: Documentation of the model.
- interactive: A script for running the model with easily-adjustable parameters and
  displaying the output.
- sensitivity\_analysis: Tools for analyzing the sensitivity of the rate profiles to the
  values for each of the input parameters.
- lake.cfg: Configuration for different tools, with lists of the names for processes and
  some default values from which to generate scripts for the various tools.
- set\_matlab\_path.sh: Convenience tool for setting the matlab path when running scripts
  from the terminal.

# Usage

You should be able to run the sensitivity analysis in parallel. For running on 20 nodes on coyote,

```
cat submit_0.sh | ssub -n 20
```