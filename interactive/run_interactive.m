% nb: be sure to run write_default_values python script first so that the
% run_interactive_defaults matlab script exists!

% include the primary matlab scripts
addpath('../bin');

run_it = 1; % run the simulation again?

show_concs = 1; % show the concentrations?
show_c = 1; % show the carbon species?

show_rates = 1; % show the rates data?

save_data = 1;  % save the time, rate, and concs data?

show_time = 100;

clf;

% grab the species map
s = species_map();

% if required, rerun the simulation
if run_it
    run_interactive_defaults;
end

if save_data
    write_data_to_files
end

if show_concs
    subplot(1,2,1);
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
    %clf;

    % show the final concentration curves
    % plot all the curves on one graph, bolding one line
    hold all;
    for i = idx
        if i <= i_bold_max
            plot(concs_history(show_time, :, i), 'LineWidth', 2)
        else
            plot(concs_history(show_time, :, i))
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
end

if show_rates
    subplot(1,2,2);
    idx = [7 8 9 10 6 5];

    %clf;
    hold all;
    for i = idx
        plot(rates_history(show_time, :, i));
    end
    hold off;

    % attach labels
    xlabel('depth');
    ylabel('rate');

    legend('aer het', 'nit red', 'iron red', 'sulf red', 'metht sulf', 'metht oxy');
end