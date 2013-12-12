% nb: be sure to run write_default_values python script first so that the
% run_interactive_defaults matlab script exists!

% include the primary matlab scripts
addpath('../bin');

run_it = 1; % run the simulation again?
show_c = 1; % show the carbon species?

% grab the species map
s = species_map();

% if required, rerun the simulation
if run_it
    run_interactive_defaults;
end

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

% show the final concentration curves
% plot all the curves on one graph, bolding one line
hold all;
for i = idx
    if i <= i_bold_max
        plot(concs_history(end, :, i), 'LineWidth', 2)
    else
        plot(concs_history(end, :, i))
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
