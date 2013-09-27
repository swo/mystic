run_it = 1;
show_o = 1;

s = species_map();

if run_it
    %[t, y, flux, bio_rates, abio_rates] = run();
    [t, c, m] = run();
end

if show_o
    to_show = {'O', 'C', 'N+', 'N-', 'S+', 'S-', 'Fe+', 'Fe-'};
else
    to_show = {'N+', 'S+', 'Fe-'};
end

idx = cellfun(@(x) s(x), to_show);

clf;

subplot(2, 1, 1)
hold all;
for i = idx
    if i == 2
        plot(c(end, :, i), 'LineWidth', 2)
    else
        plot(c(end, :, i))
    end
end
hold off;

xlabel('depth');
ylabel('concentration');

if show_o
    legend('O', 'C', 'N+', 'N-', 'S+', 'S-', 'Fe+', 'Fe-');
else
    legend('N+', 'S+', 'Fe-');
end

subplot(2, 1, 2)
hold all;
for i = 1:6
    plot(m(end, :, i))
end
hold off;

legend('resp N', 'resp Fe', 'resp S', 'ox N', 'ox S', 'ox Fe')

%flux
%bio_rates
%abio_rates