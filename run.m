function [t, y, final_flux, final_ma_op_rates, final_tea_rates] = run()

% constants
nitrogen_ratio = 0.30;  % N- released per C degraded, 0.15 from Redfield
carbon_ratio = 1.01; % C dumped per O dumped

diffusion_constant_per_compartment2 = 0.05;
fixed_oxygen_level = 50.0;
fixed_oxygen_diffusion = 1e2;

fixed_co2_level = 600;
fixed_co2_diffusion = fixed_oxygen_diffusion;

% simulation parameters
n_x = 15;
t_max = 1e4;
diffusion_constant = diffusion_constant_per_compartment2 * n_x ^ 2;
%fixed_oxygen_diffusion = diffusion_constant;

% species list
[s, ~, n_species] = species_map();

precipitation_constant_input = [
    s('N-'), 0.0
    s('S-'), 0.0
    s('Fe+'), 0.1
    s('C'), 0.1
    s('CO2'), 0.0
];

% mass action, one product reactions
ma_op_rxns = [
    % secondary oxidations
    % A + B -> C, rate constant
    % constants are from Hunter, Wang, van Cappellen 1998
    % [rates] = uM-1 yr-1
    s('O'), s('N-'), s('N+'), 5.0   % k_4^sr = 5e6 M-1 yr-1
    s('O'), s('Fe-'), s('Fe+'), 10.0    % k_2^sr = 1e7 M-1 yr-1
    s('O'), s('S-'), s('S+'), 0.16  % k_5^sr = 1.6e5 M-1 yr-1
];
[n_ma_op_rxns, ~] = size(ma_op_rxns);

% primary oxidation
% primary oxidation rate constant: HWvC report 3e-5 to 3e1 yr-1
po_rc = 1e0;

% primary oxidation terminal electron acceptors
po_teas = [
    % in species, out species, c_lim, # electrons (e_j)
    % c_lim in uM
    s('O'), s('null'), 20.0, 4  % output is water
    s('N+'), s('null'), 5.0, 5  % output is N2
    s('Fe+'), s('Fe-'), 1.0, 1 % had to adjust from HWvC on account of units (60.0)
    s('S+'), s('S-'), 0.03, 8
    s('CO2'), s('null'), 0.0, 8 % output is methane
];
[n_po_teas, ~] = size(po_teas);


% initialize the lake
concs0 = zeros(n_x, n_species);
concs0(:, s('C')) = 50.0;
concs0(:, s('O')) = 50.0;

concs0(:, s('N+')) = 100.0;
concs0(:, s('N-')) = 100.0;

concs0(:, s('Fe+')) = 30.0;
concs0(:, s('Fe-')) = 30.0;

concs0(:, s('S+')) = 60.0;
concs0(:, s('S-')) = 60.0;


% compute internal parameters

% make the precipitation constants
precipitation_constants = zeros([1 n_species]);
for i = 1: length(precipitation_constant_input)
    idx = precipitation_constant_input(i, 1);
    val = precipitation_constant_input(i, 2);

    precipitation_constants(idx) = val;
end

% make new precipitation values
D = diffusion_constant;
D_plus = (1.0 + precipitation_constants) * D;
D_minus = (1.0 - precipitation_constants) * D;

% grab the unchanging columns from the reaction matrix
ma_op_reac1_i = ma_op_rxns(:, 1);
ma_op_reac2_i = ma_op_rxns(:, 2);
ma_op_prod_i = ma_op_rxns(:, 3);
ma_op_rc = ma_op_rxns(:, 4)';

po_tea_i = po_teas(:, 1);
po_tea_prod_i = po_teas(:, 2);
po_tea_clim = po_teas(:, 3);
po_tea_e = po_teas(:, 4)';


% -- flux function --
n_total = n_x * n_species;

function [ma_op_rates, tea_rates] = rates(concs_row)
    % compute the mass action rates
    ma_op_reac1 = concs_row(ma_op_reac1_i);
    ma_op_reac2 = concs_row(ma_op_reac2_i);

    ma_op_rates = ma_op_rc .* ma_op_reac1 .* ma_op_reac2;

    % compute the fraction of electrons consumed by each TEA
    f = zeros(1, n_po_teas);
    for i = 1: n_po_teas
        f(i) = (1.0 - sum(f)) * min(1, concs_row(po_tea_i(i)) / po_tea_clim(i));
    end

    % assign whatever fraction is left to methanogenesis
    %f_methanogenesis = 1.0 - sum(f);

    % weight each TEA fraction by their # electrons
    tea_rates = po_rc * concs_row(s('C')) * f ./ po_tea_e;
end

% twiddle because this is time-independent
function [conc_fluxes] = flux(~, concs_vector)
    concs = reshape(concs_vector, [n_x, n_species]);

    conc_fluxes = zeros(n_x, n_species);
    
    % apply the fixed oxygen term
    oxygen_source = fixed_oxygen_diffusion * (fixed_oxygen_level - concs(1, s('O')));
    conc_fluxes(1, s('O')) = conc_fluxes(1, s('O')) + oxygen_source;
    conc_fluxes(1, s('C')) = conc_fluxes(1, s('C')) + carbon_ratio * oxygen_source;
    
    % apply the fixed co2 term
    co2_source = fixed_co2_diffusion * (fixed_co2_level - concs(1, s('CO2')));
    conc_fluxes(1, s('CO2')) = conc_fluxes(1, s('CO2')) + co2_source;
    
    for x = 1: n_x
        [ma_op_rates, tea_rates] = rates(concs(x, :));

        % apply the mass action rates
        conc_fluxes(x, :) = conc_fluxes(x, :) - accumarray(ma_op_reac1_i, ma_op_rates, [n_species, 1])';
        conc_fluxes(x, :) = conc_fluxes(x, :) - accumarray(ma_op_reac2_i, ma_op_rates, [n_species, 1])';
        conc_fluxes(x, :) = conc_fluxes(x, :) + accumarray(ma_op_prod_i, ma_op_rates, [n_species, 1])';

        % apply the primary oxidation rates
        conc_fluxes(x, :) = conc_fluxes(x, :) - accumarray(po_tea_i, tea_rates, [n_species, 1])';
        conc_fluxes(x, :) = conc_fluxes(x, :) + accumarray(po_tea_prod_i, tea_rates, [n_species, 1])';
        
        total_tea = sum(tea_rates);
        conc_fluxes(x, s('C')) = conc_fluxes(x, s('C')) - total_tea;
        conc_fluxes(x, s('CO2')) = conc_fluxes(x, s('CO2')) + total_tea;
        conc_fluxes(x, s('N-')) = conc_fluxes(x, s('N-')) + nitrogen_ratio * total_tea;
        
        % -- diffusion --        
        if x > 1
            conc_fluxes(x, :) = conc_fluxes(x, :) + D_plus .* concs(x - 1, :) - D_minus .* concs(x, :);
        end

        if x < n_x
            conc_fluxes(x, :) = conc_fluxes(x, :) - D_plus .* concs(x, :) + D_minus .* concs(x + 1, :);
        end

    end % for x
    
    conc_fluxes(:, s('null')) = 0.0;
    conc_fluxes = reshape(conc_fluxes, [n_total, 1]);
end


% -- run the ode solver --
options = odeset('NonNegative', 1: n_total);
concs0_vector = reshape(concs0, [n_total, 1]);
[t, y] = ode15s(@flux, [0 t_max], concs0_vector, options);

[n_time_slices, ~] = size(y);
y = reshape(y, n_time_slices, n_x, n_species);

concs_vector = reshape(y(end, :, :), [n_total, 1]);
final_flux = reshape(flux(0, concs_vector), [n_x, n_species]);

final_ma_op_rates = zeros(n_x, n_ma_op_rxns);
final_tea_rates = zeros(n_x, n_po_teas);
for x = 1: n_x
    [final_ma_op_rates(x, :), final_tea_rates(x, :)] = rates(squeeze(y(end, x, :))');
end

end
