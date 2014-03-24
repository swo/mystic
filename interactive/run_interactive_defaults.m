[time_slices, concs_history, rates_history] = lake( ...
	0.0,	... oxygen_bubble_rate
	0.0,	... nitrogen_source
	0.1,	... nitrogen_ratio
	4.5e4,	... carbon_source
	0.5e4,	... oxygen_source
	1500.0,	... methane_source
	1.0,	... t_max
	0.3,	... fe_precipitation
    0.1,    ... carbon precip
	50.0,	... diffusion_constant
	1.0e4,	... ma_op_o_fe_rate_const
	5.0,	... ma_op_o_n_rate_const
	0.16,	... ma_op_o_s_rate_const
	0.01,	... ma_op_fe_n_rate_const
	1000.0,	... ma_op_ch4_o_rate_const
	0.01,	... ma_op_ch4_s_rate_const
	1,	... primary_ox_rate_const
	20.0,	... c_lim_o
	5.0,	... c_lim_n
	0.1,	... c_lim_fe
	30,	... c_lim_s
	200.0,	... concs0_c
	50.0,	... concs0_o
	100.0,	... concs0_ntot
	1.0,	... pm_ratio_n
	60.0,	... concs0_fetot
	1.0,	... pm_ratio_fe
	200.0,	... concs0_stot
	1.0	... pm_ratio_s
);