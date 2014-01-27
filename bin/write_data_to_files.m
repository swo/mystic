% assume that rates_history, time_slices, and concs_history are loaded already

prefix = 'data';

% save the time data
fn = sprintf('%s/time.txt', prefix);
dlmwrite(fn, time_slices);

% save the rate data
[~, ~, n_rates] = size(rates_history);
for i = 1: n_rates
    fn = sprintf('%s/rate_%01d.txt', prefix, i);
    dlmwrite(fn, squeeze(rates_history(:, :, i))', 'delimiter', '\t');
end

% save the concentration profiles
[~, ~, n_concs] = size(concs_history);
for i = 1: n_concs
    fn = sprintf('%s/concs_%01d.txt', prefix, i);
    dlmwrite(fn, squeeze(concs_history(:, :, i))', 'delimiter', '\t');
end

% save all of them together for good measure
fn = sprintf('%s/history.mat', prefix);
save(fn, 'time_slices', 'rates_history', 'concs_history');
