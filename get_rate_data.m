% pick evenly spaced time points
n_t = 50;
t_min = 1;
t_max = max(t);
ts = linspace(t_min, t_max, n_t);

clearvars('idx');
for i = 1: n_t
    this_time = ts(i);
    [~, this_index] = min(abs(t - this_time));
    idx(i) = this_index;
end
idx

% save the time data
dlmwrite('rate_data/time.csv', t(idx));

% write the mass action data
[~, ~, n_r] = size(final_ma_op_rates);
for r = 1: n_r
    fn = sprintf('rate_data/ma_%01d.csv', r);
    dlmwrite(fn, squeeze(final_ma_op_rates(idx, :, r))');
end

% write the TEA data
[~, ~, n_r] = size(final_tea_rates);
for r = 1: n_r
    fn = sprintf('rate_data/tea_%01d.csv', r);
    dlmwrite(fn, squeeze(final_tea_rates(idx, :, r))');
end

% write the concentration profiles
[~, ~, n_c] = size(y);
for c = 1: n_c
    fn = sprintf('rate_data/conc_%01d.csv', c);
    dlmwrite(fn, squeeze(y(idx, :, c))');
end

% write out the names of the species in the concentration profiles
[~, species, n_s] = species_map();
fh = fopen('rate_data/conc_names.csv', 'wt');
for i = 1: n_s
    fprintf(fh, '%s\n', species{i});
end
fclose(fh);