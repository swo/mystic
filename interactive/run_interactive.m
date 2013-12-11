% include the primary matlab scripts
addpath('../bin');

run_it = 1; % run the simulation again?
show_c = 1; % show the carbon species?

% grab the species map
s = species_map();

% if required, rerun the simulation
if run_it
    [t, y, final_flux, final_ma_op_rates, final_tea_rates] = run();
end

% spit out the final results from the model for cursory diagnostics
final_flux
final_ma_op_rates
final_tea_rates

% prepare lists of the curves to be drawn
if show_c
    to_show = {'O', 'N+', 'N-', 'S+', 'S-', 'Fe+', 'Fe-', 'CH4'};
    i_bold_max = 2;
else
    to_show = {'O', 'N+', 'N-', 'S+', 'S-', 'Fe+', 'Fe-'};
    i_bold_max = 1;
end
idx = cellfun(@(x) s(x), to_show);

% clear the plot
clf;

% plot all the curves on one graph, bolding one line
hold all;
for i = idx
    if i <= i_bold_max
        plot(y(end, :, i), 'LineWidth', 2)
    else
        plot(y(end, :, i))
    end
end
hold off;

% attach labels
xlabel('depth');
ylabel('concentration');

% attach legend
if show_c
    legend('O', 'N+', 'N-', 'S+', 'S-', 'Fe+', 'Fe-', 'CH4');
else
    legend('O', 'N+', 'N-', 'S+', 'S-', 'Fe+', 'Fe-');
end
