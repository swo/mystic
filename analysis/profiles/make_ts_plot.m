function [] = make_ts_plot(dat, fn, n_exclude, t_max)
%MAKE_TS_PLOT Make a timeseries heatmap plot
%   Takes care of the font sizes, etc.
%   n_exclude : int
%       number of timepoints at the start to skip
%   t_max : float
%       scale for the time axis


[dat_max, ~] = size(dat);

imagesc(dat(n_exclude + 1:end, :)');
cb = colorbar;

set(gca, 'FontSize', 30);
set(cb, 'FontSize', 30);

L = get(gca, 'XLim');
ticks = linspace(L(1), L(2), 4);
set(gca, 'XTick', ticks);
set(gca, 'XTickLabel', floor((ticks + n_exclude) * t_max / dat_max));
set(gca, 'YTick', [0.5, 5, 10, 15]);
set(gca, 'YTickLabel', {'5', '10', '15', '20'});

L = get(cb, 'YLim');
ticks = linspace(L(1), L(2), 3);
set(cb, 'YTick', ticks);
set(cb, 'YTickLabel', max(0, floor(ticks)));

print(fn, '-dpdf');

end
