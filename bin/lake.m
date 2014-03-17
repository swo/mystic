function [time_slices, concs_history, rates_history] = lake(NITROGEN_RATIO, FIXED_CARBON_LEVEL, FIXED_OXYGEN_LEVEL, FIXED_OXYGEN_DIFFUSION, FIXED_BOTTOM_METHANE, T_MAX, FE_PRECIPITATION, DIFF_CONST_COMP, MA_OP_O_FE_RATE_CONST, MA_OP_O_N_RATE_CONST, MA_OP_O_S_RATE_CONST, MA_OP_FE_N_RATE_CONST, MA_OP_CH4_O_RATE_CONST, MA_OP_CH4_S_RATE_CONST, PRIMARY_OX_RATE_CONST, C_LIM_O, C_LIM_N, C_LIM_FE, C_LIM_S, CONCS0_C, CONCS0_O, CONCS0_NTOT, PM_RATIO_N, CONCS0_FETOT, PM_RATIO_FE, CONCS0_STOT, PM_RATIO_S)
%% Constants
% These are constants that make assertions about the actual system

nitrogen_ratio = NITROGEN_RATIO;  % N- released per C degraded, 0.15 from Redfield

diffusion_constant_per_compartment2 = DIFF_CONST_COMP; % input diffusion constant
fixed_oxygen_level = FIXED_OXYGEN_LEVEL;  % oxygen level at thermocline
fixed_oxygen_diffusion = FIXED_OXYGEN_DIFFUSION;   % diffusion from oxygen above the thermocline
fixed_oxygen_compartments = 1;

fixed_top_methane_level = 0.0;
fixed_bottom_methane_level = FIXED_BOTTOM_METHANE;
fixed_methane_diffusion = fixed_oxygen_diffusion;

%% Simulation parameters
% These are constants that affect the simulation but do not make
% assertions about the actual system.

n_x = 17;   % number of compartments
t_max = T_MAX;   % time until end of simulation (yrs)
n_time_slices = 100;

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
    s('Fe+'), FE_PRECIPITATION
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
    0.25, s('O'), 1, s('Fe-'), 1, s('Fe+'), MA_OP_O_FE_RATE_CONST    % k_2^sr = 1e7 M-1 yr-1
    2, s('O'), 1, s('N-'), 1, s('N+'), MA_OP_O_N_RATE_CONST   % k_4^sr = 5e6 M-1 yr-1
    2, s('O'), 1, s('S-'), 1, s('S+'), MA_OP_O_S_RATE_CONST  % k_5^sr = 1.6e5 M-1 yr-1
    5, s('Fe-'), 1, s('N+'), 5, s('Fe+'), MA_OP_FE_N_RATE_CONST   % swo> my guess
    1, s('CH4'), 2, s('O'), 1, s('null'), MA_OP_CH4_O_RATE_CONST    % k_9^sr = 1e10 M-1 yr-1
    1, s('CH4'), 1, s('S+'), 1, s('S-'), MA_OP_CH4_S_RATE_CONST    % k_10^sr = 1e5 M-1 yr-1
];
[n_ma_op_rxns, ~] = size(ma_op_rxns);

% primary oxidation reactions
% primary oxidation rate constant (po_rc): HWvC report 3e-5 to 3e1 yr-1
po_rc = PRIMARY_OX_RATE_CONST;

% primary oxidation terminal electron acceptors (po_teas)
po_teas = [
    % columns: in species, out species, c_lim, # electrons (e_j)
    % units: [c_lim] = uM
    s('O'), s('null'), C_LIM_O, 4  % HWcV O2_lim=20; output is water
    s('N+'), s('null'), C_LIM_N, 5  % 5; output is N2
    s('Fe+'), s('Fe-'), C_LIM_FE, 1 % 0.1; had to adjust from HWvC on account of units (60.0)
    s('S+'), s('S-'), C_LIM_S, 8 % note HWvC have 0.03 mM (= 30 uM)
    s('null'), s('CH4'), 0.0, 8 % output is methane
];
[n_po_teas, ~] = size(po_teas);

%% Initial concentrations
% initialize the lake, asserting flat profiles for each metabolite
% metabolites not mentioned have concentration 0
concs0 = zeros(n_x, n_species);

concs0(:, s('C')) = CONCS0_C;
concs0(:, s('O')) = CONCS0_O;

%Start with an initial conc of N (CONCS0_NTOT) and a ratio of plus to minus PM_RATIO_N
CONCS0_NPLUS = CONCS0_NTOT/(1 + (1/PM_RATIO_N));
CONCS0_NMINUS = CONCS0_NTOT/(PM_RATIO_N + 1);
concs0(:, s('N+')) = CONCS0_NPLUS;
concs0(:, s('N-')) = CONCS0_NMINUS;

%Start with an initial conc of FE (CONCS0_FETOT) and a ratio of plus to minus PM_RATIO_FE
CONCS0_FEPLUS = CONCS0_FETOT/(1 + (1/PM_RATIO_FE));
CONCS0_FEMINUS = CONCS0_FETOT/(PM_RATIO_FE + 1);
concs0(:, s('Fe+')) = CONCS0_FEPLUS;
concs0(:, s('Fe-')) = CONCS0_FEMINUS;

%Start with an initial conc of S (CONCS0_STOT) and a ratio of plus to minus PM_RATIO_S
CONCS0_SPLUS = CONCS0_STOT/(1 + (1/PM_RATIO_S));
CONCS0_SMINUS = CONCS0_STOT/(PM_RATIO_S + 1);
concs0(:, s('S+')) = CONCS0_SPLUS;
concs0(:, s('S-')) = CONCS0_SMINUS;


%% Internal parameters
% These are parameters that simplify the code but make no assertions about
% the mechanics of the simulation or the actual system

% make the precipitation constants
precipitation_constants = zeros([1 n_species]);
for i = 1: size(precipitation_constant_input, 1)
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

function [ma_op_rates, po_carbon_rate, tea_rates] = rates(concs_row)
    % compute the mass action rates
    ma_op_reac1 = concs_row(ma_op_reac1_i);
    ma_op_reac2 = concs_row(ma_op_reac2_i);

    ma_op_rates = ma_op_rc .* ma_op_reac1 .* ma_op_reac2;

    % compute the fraction of electrons consumed by each TEA, then give the rest to methanogenesis
    f = zeros(1, n_po_teas);
    for i = 1: n_po_teas - 1
        f(i) = (1.0 - sum(f)) * min(1, concs_row(po_tea_i(i)) / po_tea_clim(i));
    end
    f(end) = 1.0 - sum(f);

    % weight each TEA fraction by their # electrons
    po_carbon_rate = po_rc * concs_row(s('C'));
    tea_rates = po_carbon_rate * f ./ po_tea_e;
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
    oxygen_source = fixed_oxygen_diffusion * (fixed_oxygen_level - concs(fixed_oxygen_compartments, s('O')));
    conc_fluxes(fixed_oxygen_compartments, s('O')) = conc_fluxes(fixed_oxygen_compartments, s('O')) + oxygen_source;

    % apply fixed carbon
    % swo> used oxygen here, was lazy
    carbon_source = fixed_oxygen_diffusion * (FIXED_CARBON_LEVEL - concs(fixed_oxygen_compartments, s('C')));
    conc_fluxes(fixed_oxygen_compartments, s('C')) = conc_fluxes(fixed_oxygen_compartments, s('C')) + carbon_source;
    
    % apply the fixed methane level at the thermocline
    methane_source = fixed_methane_diffusion * (fixed_top_methane_level - concs(1, s('CH4')));
    conc_fluxes(1, s('CH4')) = conc_fluxes(1, s('CH4')) + methane_source;

    % apply the methane source at the bottom of the lake
    methane_source = fixed_methane_diffusion * (fixed_bottom_methane_level - concs(n_x, s('CH4')));
    conc_fluxes(n_x, s('CH4')) = conc_fluxes(n_x, s('CH4')) + methane_source;
    
    for x = 1: n_x
        [ma_op_rates, po_carbon_rate, tea_rates] = rates(concs(x, :));

        % apply the mass action rates
        conc_fluxes(x, :) = conc_fluxes(x, :) - accumarray(ma_op_reac1_i, ma_op_reac1_c .* ma_op_rates, [n_species, 1])';
        conc_fluxes(x, :) = conc_fluxes(x, :) - accumarray(ma_op_reac2_i, ma_op_reac2_c .* ma_op_rates, [n_species, 1])';
        conc_fluxes(x, :) = conc_fluxes(x, :) + accumarray(ma_op_prod_i, ma_op_prod_c .* ma_op_rates, [n_species, 1])';

        % apply the primary oxidation rates
        conc_fluxes(x, :) = conc_fluxes(x, :) - accumarray(po_tea_i, tea_rates, [n_species, 1])';
        conc_fluxes(x, :) = conc_fluxes(x, :) + accumarray(po_tea_prod_i, tea_rates, [n_species, 1])';
        
        conc_fluxes(x, s('C')) = conc_fluxes(x, s('C')) - po_carbon_rate;
        conc_fluxes(x, s('N-')) = conc_fluxes(x, s('N-')) + nitrogen_ratio * po_carbon_rate;
        
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
[time_slices, y] = ode15s(@flux, linspace(0.0, t_max, n_time_slices), concs0_vector, options);

% unfold the result y, putting it into a 3D space whose dimensions
% correspond to time, depth, and metabolite
[n_time_slices, ~] = size(y);
concs_history = reshape(y, n_time_slices, n_x, n_species);

%% get all the reaction rates for all timepoints
% ma_op then teas
rates_history = zeros(n_time_slices, n_x, n_ma_op_rxns + n_po_teas);
for time = 1: n_time_slices
    for x = 1: n_x
        [rates_history(time, x, 1:n_ma_op_rxns), ~, rates_history(time, x, n_ma_op_rxns + 1:end)] = rates(squeeze(concs_history(time, x, :))');
    end
end

end
