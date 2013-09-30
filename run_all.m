run_it = 1;
show_o = 1;

s = species_map();

if run_it
    %[t, y, flux, bio_rates, abio_rates] = run();
    [t, c, m] = run();
end

if show_o
    to_show = {'O', 'C', 'N+', 'N-', 'S+', 'S-', 'Fe+', 'Fe-', 'P'};
else
    to_show = {'N+', 'S+', 'Fe-'};
end

idx = cellfun(@(x) s(x), to_show);

clf;

subplot(2, 1, 1)
hold all;
for i = idx
    if i == 2 || i == 3
        plot(c(end, :, i), 'LineWidth', 2)
    else
        plot(c(end, :, i))
    end
end
hold off;

xlabel('depth');
ylabel('concentration');

if show_o
    legend('O', 'C', 'N+', 'N-', 'S+', 'S-', 'Fe+', 'Fe-', 'P');
else
    legend('N+', 'S+', 'Fe-');
end

subplot(2, 1, 2)
hold all;
for i = 1:7
    plot(m(end, :, i))
end
hold off;

legend('aer het', 'denit', 'resp Fe', 'resp S', 'ox N', 'ox S', 'NFE')

%flux
%bio_rates
%abio_rates