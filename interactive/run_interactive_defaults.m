[time_slices, concs_history, rates_history] = run(
	0.3,	... nitrogen_ratio
	1.0,	... carbon_ratio
	75.0,	... fixed_oxygen_level
	10000.0,	... fixed_oxygen_diffusion
	600.0,	... fixed_co2_level
	25.0,	... t_max
	0.1,	... fe_precipitation
	0.75,	... diff_const_comp
	10.0,	... ma_op_o_fe_rate_const
	5.0,	... ma_op_o_n_rate_const
	0.16,	... ma_op_o_s_rate_const
	0.01,	... ma_op_fe_n_rate_const
	1,	... primary_ox_rate_const
	20.0,	... c_lim_o
	5.0,	... c_lim_n
	0.1,	... c_lim_fe
	30,	... c_lim_s
	0.0,	... c_lim_co2
	200.0,	... concs0_c
	50.0,	... concs0_o
	200.0,	... concs0_ntot
	1.0,	... pm_ratio_n
	60.0,	... concs0_fetot
	1.0,	... pm_ratio_fe
	240.0,	... concs0_stot
	1.0	... pm_ratio_s
);