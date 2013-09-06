run_it = 1;
show_o = 1;

s = species_map();

if run_it
    [t, y, flux, bio_rates, abio_rates] = run();
end

if show_o
    to_show = {'O', 'C', 'N+', 'N-', 'S+', 'S-', 'Fe+', 'Fe-'};
else
    to_show = {'N+', 'N-', 'S+', 'S-', 'Fe+', 'Fe-'};
end

idx = cellfun(@(x) s(x), to_show);

clf;

twod = 1;

hold all;
for i = idx
    plot(y(end, :, i))
end

xlabel('depth');
ylabel('time');
zlabel('concentration');

if show_o
    legend('O', 'C', 'N+', 'N-', 'S+', 'S-', 'Fe+', 'Fe-');
else
    legend('N+', 'N-', 'S+', 'S-', 'Fe+', 'Fe-');
end

hold off;

flux
bio_rates
abio_rates