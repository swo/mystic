function [] = run(NITROGEN_RATIO, CARBON_RATIO, FIXED_OXYGEN_LEVEL, FIXED_OXYGEN_DIFFUSION, FIXED_CO2_LEVEL, T_MAX, FE_PRECIPITATION, OUT)

%% Constants
% These are constants that make assertions about the actual system

nitrogen_ratio = NITROGEN_RATIO;  % N- released per C degraded, 0.15 from Redfield
carbon_ratio = CARBON_RATIO; % C dumped per O dumped

diffusion_constant_per_compartment2 = 0.75; % input diffusion constant
fixed_oxygen_level = FIXED_OXYGEN_LEVEL;  % oxygen level at thermocline
fixed_oxygen_diffusion = FIXED_OXYGEN_DIFFUSION;   % diffusion from oxygen above the thermocline

fixed_co2_level = FIXED_CO2_LEVEL;  % CO2 level at thermocline
fixed_co2_diffusion = fixed_oxygen_diffusion;

fixed_methane_level = 0.0;
fixed_methane_diffusion = fixed_oxygen_diffusion;

%% Simulation parameters
% These are constants that affect the simulation but do not make
% assertions about the actual system.

n_x = 17;   % number of compartments
t_max = T_MAX;   % time until end of simulation (yrs)

% compute the diffusion constant, which is dependent on the length scale,
% which depends on the square of the number of compartments
diffusion_constant = diffusion_constant_per_compartment2 * n_x ^ 2;

%% Species map
% import the species list using the separate function file:
% "s" is a hash from a string that names the species to its index in the
% concentration matrix
[s, ~, n_species] = species_map();

%% Precipitation constants
% assert the precipitation constants: a constant 0 means equal diffusion up
% and down; +0.1 means a molecule is 10% as likely to go down as to go up;
% -0.1 means 10% more likley to go up than down
precipitation_constant_input = [
    s('N-'), 0.0
    s('S-'), 0.0
    s('Fe+'), FE_PRECIPITATION
    s('C'), 0.0
    s('CO2'), 0.0
];

%% Reaction constants
% there are mass action reactions and the primary oxidations
% primary oxidations follow Hunter, Wang, and Van Cappellen 1998
% (abbreviated HWvC)

% mass action, one product reactions (ma_op_rxns)
ma_op_rxns = [
    % secondary oxidations
    % aA + bB -> cC (columns: a, A, b, B, c, C, rate constant)
    % constants k_i are from HWvC
    % units: [rates] = uM-1 yr-1
    0.25, s('O'), 1, s('Fe-'), 1, s('Fe+'), 10.0    % k_2^sr = 1e7 M-1 yr-1
    2, s('O'), 1, s('N-'), 1, s('N+'), 5.0   % k_4^sr = 5e6 M-1 yr-1
    2, s('O'), 1, s('S-'), 1, s('S+'), 0.16  % k_5^sr = 1.6e5 M-1 yr-1
    5, s('Fe-'), 1, s('N+'), 5, s('Fe+'), 0.01   % swo> my guess
];
[n_ma_op_rxns, ~] = size(ma_op_rxns);

% primary oxidation reactions
% primary oxidation rate constant (po_rc): HWvC report 3e-5 to 3e1 yr-1
po_rc = 1e0;

% primary oxidation terminal electron acceptors (po_teas)
po_teas = [
    % columns: in species, out species, c_lim, # electrons (e_j)
    % units: [c_lim] = uM
    s('O'), s('null'), 20.0, 2  % HWcV O2_lim=20; output is water
    s('N+'), s('null'), 5.0, 5  % 5; output is N2
    s('Fe+'), s('Fe-'), 0.1, 1 % 0.1; had to adjust from HWvC on account of units (60.0)
    s('S+'), s('S-'), 30, 8 % note HWvC have 0.03 mM (= 30 uM)
    s('CO2'), s('CH4'), 0.0, 8 % output is methane
];
[n_po_teas, ~] = size(po_teas);

%% Initial concentrations
% initialize the lake, asserting flat profiles for each metabolite
% metabolites not mentioned have concentration 0
concs0 = zeros(n_x, n_species);
concs0(:, s('C')) = 200.0;
concs0(:, s('O')) = 50.0;

concs0(:, s('N+')) = 100.0;
concs0(:, s('N-')) = 100.0;

concs0(:, s('Fe+')) = 30.0;
concs0(:, s('Fe-')) = 30.0;

concs0(:, s('S+')) = 120.0;
concs0(:, s('S-')) = 120.0;


%% Internal parameters
% These are parameters that simplify the code but make no assertions about
% the mechanics of the simulation or the actual system

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
ma_op_reac1_c = ma_op_rxns(:, 1)';
ma_op_reac1_i = ma_op_rxns(:, 2);
ma_op_reac2_c = ma_op_rxns(:, 3)';
ma_op_reac2_i = ma_op_rxns(:, 4);
ma_op_prod_c = ma_op_rxns(:, 5)';
ma_op_prod_i = ma_op_rxns(:, 6);
ma_op_rc = ma_op_rxns(:, 7)';

po_tea_i = po_teas(:, 1);
po_tea_prod_i = po_teas(:, 2);
po_tea_clim = po_teas(:, 3);
po_tea_e = po_teas(:, 4)';


%% Define the flux functions

% -- rates --
% This function takes a row from the concentration matrix (i.e., a
% horizontal slice from the lake) and computes the rates of the mass action
% and primary oxidation reactions

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

% -- flux --
% This function is taken as an argument by the ODE solver, which feeds this
% function the concentration matrix (flattened into a concentration
% vector). The function computes the rates of reactions, diffusions, and
% precipitations and outputs the fluxes for each metabolite at each depth.

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
    
    methane_source = fixed_methane_diffusion * (fixed_methane_level - concs(1, s('CH4')));
    conc_fluxes(1, s('CH4')) = conc_fluxes(1, s('CH4')) + methane_source;
    
    for x = 1: n_x
        [ma_op_rates, tea_rates] = rates(concs(x, :));

        % apply the mass action rates
        conc_fluxes(x, :) = conc_fluxes(x, :) - accumarray(ma_op_reac1_i, ma_op_reac1_c .* ma_op_rates, [n_species, 1])';
        conc_fluxes(x, :) = conc_fluxes(x, :) - accumarray(ma_op_reac2_i, ma_op_reac2_c .* ma_op_rates, [n_species, 1])';
        conc_fluxes(x, :) = conc_fluxes(x, :) + accumarray(ma_op_prod_i, ma_op_prod_c .* ma_op_rates, [n_species, 1])';

        % apply the primary oxidation rates
        conc_fluxes(x, :) = conc_fluxes(x, :) - accumarray(po_tea_i, tea_rates, [n_species, 1])';
        conc_fluxes(x, :) = conc_fluxes(x, :) + accumarray(po_tea_prod_i, tea_rates, [n_species, 1])';
        
        total_tea = sum(tea_rates);
        conc_fluxes(x, s('C')) = conc_fluxes(x, s('C')) - total_tea;
        conc_fluxes(x, s('CO2')) = conc_fluxes(x, s('CO2')) + total_tea;
        conc_fluxes(x, s('N-')) = conc_fluxes(x, s('N-')) + nitrogen_ratio * total_tea;
        
        % diffusion      
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


%% ODE solver
% This section feeds the flux function to the ODE solver.

% all concentrations are constrained to be nonnegative
options = odeset('NonNegative', 1: n_total);

% initially flatten the concentration matrix
concs0_vector = reshape(concs0, [n_total, 1]);

% run the ODE solver (ode15s)
% t is the times at which the ODE solver gives output. They are not evenly
% spaced! y is a matrix whose rows are the flattened concentration matrices
% at each time step
[t, y] = ode15s(@flux, [0 t_max], concs0_vector, options);

% unfold the result y, putting it into a 3D space whose dimensions
% correspond to time, depth, and metabolite
[n_time_slices, ~] = size(y);
y = reshape(y, n_time_slices, n_x, n_species);

% get the concentration profile at the last timepoint, flatten it to a
% vector, and feed that to the flux function to get the fluxes at the last
% timepoint
concs_vector = reshape(y(end, :, :), [n_total, 1]);
final_flux = reshape(flux(0, concs_vector), [n_x, n_species]);

% similarly, get the reaction rates at the last timepoint
final_ma_op_rates = zeros(n_x, n_ma_op_rxns);
final_tea_rates = zeros(n_x, n_po_teas);
for x = 1: n_x
    [final_ma_op_rates(x, :), final_tea_rates(x, :)] = rates(squeeze(y(end, x, :))');
end

all_rates = [final_ma_op_rates final_tea_rates];
csvwrite(OUT, all_rates);

end
