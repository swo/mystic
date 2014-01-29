% load up the data
load('history.mat');

[n_t, n_d, n_r] = size(rates_history);
divergences = zeros(n_d, n_d);

for i = 1: n_d
    for j = 1: n_d
        % get the rows
        r1 = rates_history(end, i, :);
        r2 = rates_history(end, j, :);

        % normalize rows
        r1 = r1 / sum(r1);
        r2 = r2 / sum(r2);

        % get arthmetic mean
        avg = 0.5 * (r1 + r2);

        % removes zeros
        r1 = r1(r1 > 0.0);
        r2 = r2(r2 > 0.0);
        avg = avg(avg > 0.0);

        % get the divergence, equally weighted
        divergences(i, j) = -sum(avg .* log(avg)) - 0.5 * (-sum(r1 .* log(r1)) + -sum(r2 .* log(r2)));
    end
end

h = pcolor(divergences);
axis square
set(h, 'edgecolor', 'none');

print('divergences.pdf', '-dpdf')
