run_it = 1;
show_c = 1;

s = species_map();

if run_it
    [t, y, final_flux, final_ma_op_rates, final_tea_rates] = run();
end

final_flux
final_ma_op_rates
final_tea_rates

if show_c
    to_show = {'O', 'C', 'N+', 'N-', 'S+', 'S-', 'Fe+', 'Fe-'};
    i_bold_max = 2;
else
    to_show = {'O', 'N+', 'N-', 'S+', 'S-', 'Fe+', 'Fe-', 'CO2'};
    i_bold_max = 1;
end

idx = cellfun(@(x) s(x), to_show);

clf;

%subplot(2, 1, 1)
hold all;
for i = idx
    if i <= i_bold_max
        plot(y(end, :, i), 'LineWidth', 2)
    else
        plot(y(end, :, i))
    end
end
hold off;

xlabel('depth');
ylabel('concentration');

if show_c
    legend('O', 'C', 'N+', 'N-', 'S+', 'S-', 'Fe+', 'Fe-', 'CO2');
else
    legend('O', 'N+', 'N-', 'S+', 'S-', 'Fe+', 'Fe-', 'CO2');
end

%subplot(2, 1, 2)
%hold all;
%for i = 1:7
%    plot(m(end, :, i))
%end
%hold off;

%legend('aer het', 'denit', 'resp Fe', 'resp S', 'ox N', 'ox S', 'NFE')

%flux
%bio_rates
%abio_rates