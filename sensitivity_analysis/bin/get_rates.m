function [] = get_rates(nitrogen_ratio, carbon_ratio, fixed_oxygen_level, fixed_oxygen_diffusion, fixed_co2_level, fixed_top_methane, fixed_bottom_methane, t_max, fe_precipitation, diff_const_comp, ma_op_o_fe_rate_const, ma_op_o_n_rate_const, ma_op_o_s_rate_const, ma_op_fe_n_rate_const, ma_op_ch4_o_rate_const, ma_op_ch4_s_rate_const, primary_ox_rate_const, c_lim_o, c_lim_n, c_lim_fe, c_lim_s, c_lim_co2, concs0_c, concs0_o, concs0_ntot, pm_ratio_n, concs0_fetot, pm_ratio_fe, concs0_stot, pm_ratio_s, out_fn)

[time_slices, concs_history, rates_history] = run(nitrogen_ratio, carbon_ratio, fixed_oxygen_level, fixed_oxygen_diffusion, fixed_co2_level, fixed_top_methane, fixed_bottom_methane, t_max, fe_precipitation, diff_const_comp, ma_op_o_fe_rate_const, ma_op_o_n_rate_const, ma_op_o_s_rate_const, ma_op_fe_n_rate_const, ma_op_ch4_o_rate_const, ma_op_ch4_s_rate_const, primary_ox_rate_const, c_lim_o, c_lim_n, c_lim_fe, c_lim_s, c_lim_co2, concs0_c, concs0_o, concs0_ntot, pm_ratio_n, concs0_fetot, pm_ratio_fe, concs0_stot, pm_ratio_s);
final_rates = rates_history(end, :, :);
dlmwrite(out_fn, final_rates);

end