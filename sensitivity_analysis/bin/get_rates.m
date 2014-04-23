function [] = get_rates(oxygen_bubble_rate, nitrogen_source, nitrogen_ratio, carbon_source, oxygen_source, methane_source, t_max, fe_precipitation, carbon_precipitation, diffusion_constant, ma_op_o_fe_rate_const, ma_op_o_n_rate_const, ma_op_o_s_rate_const, ma_op_fe_n_rate_const, ma_op_ch4_o_rate_const, ma_op_ch4_s_rate_const, primary_ox_rate_const, c_lim_o, c_lim_n, c_lim_fe, c_lim_s, concs0_c, concs0_o, concs0_ntot, pm_ratio_n, concs0_fetot, pm_ratio_fe, concs0_stot, pm_ratio_s, rates_out_fn, concs_out_fn)

[time_slices, concs_history, rates_history] = lake(oxygen_bubble_rate, nitrogen_source, nitrogen_ratio, carbon_source, oxygen_source, methane_source, t_max, fe_precipitation, carbon_precipitation, diffusion_constant, ma_op_o_fe_rate_const, ma_op_o_n_rate_const, ma_op_o_s_rate_const, ma_op_fe_n_rate_const, ma_op_ch4_o_rate_const, ma_op_ch4_s_rate_const, primary_ox_rate_const, c_lim_o, c_lim_n, c_lim_fe, c_lim_s, concs0_c, concs0_o, concs0_ntot, pm_ratio_n, concs0_fetot, pm_ratio_fe, concs0_stot, pm_ratio_s);

final_rates = squeeze(rates_history(end, :, :));
dlmwrite(rates_out_fn, final_rates, precision=15);

final_concs = squeeze(concs_history(end, :, :));
dlmwrite(concs_out_fn, final_concs, precision=15);

end