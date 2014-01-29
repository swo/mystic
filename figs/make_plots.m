% load up the data
load('history.mat')
conc_names = textread('conc_names.txt', '%s', 'delimiter', '\n');
rate_names = textread('rate_names.txt', '%s', 'delimiter', '\n');
depths = dlmread('depths.txt');

% figure out the number of figures to make
[~, ~, n_concs] = size(concs_history);
[~, ~, n_rates] = size(rates_history);

% make each of the concs figures
for i = 1: n_concs
    fn = sprintf('plots/end/c_%s.pdf', conc_names{i});
    plot(depths, concs_history(end, :, i));
    print(fn, '-dpdf');

    fn = sprintf('plots/time/c_%s.pdf', concs_names{i});
    imagesc(concs_history(:, :, i));
    print(fn, '-dpdf');
end

% make each of the rates figures
for i = 1: n_rates
    fn = sprintf('plots/end/r_%s.pdf', rate_names{i});
    plot(depths, rates_history(end, :, i));
    print(fn, '-dpdf');

    fn = sprintf('plots/time/r_%s.pdf', rates_names{i});
    imagesc(rates_history(:, :, i));
    print(fn, '-dpdf');
end
