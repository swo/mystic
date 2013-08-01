% todo:
% add an oxygen production term: oxygen can only be produced near the
% surface, but it requires carbon species that flow from the top

% think again about the carbon cycle

function [sol] = run()

% constants
diffusion_constant = 0.1;
RT = 2.49e-6; % 8.3 J K^-1 mol^-1 * 300 K = 2.49 kJ mol^-1 = 2.49e-6 kJ umol^-1
rate_constant = 1e4;
faraday_constant = 9.6485e-5; % 96.4 kJ volt^-1 mol^-1 = 96.4e-6 umol^-1

% assertive parameters
photon_depth_scale = 1.0; % 1/e distance for photosynthesis (meters)
photo_delta_G_standard = -1e-4;

% methanogenesis parameters
mg_rate_constant = 1.0;
mg_delta_G_modifier = 0.0; % accounts for the [H20]^2/[H2]^4 in Q

metabolic_cutoff = 0; % Canfield's cutoff for useful metabolism; -20 kJ mol^-1 = -2e-5 kJ mmol^-1

% simulation parameters
x_max = 15;
x_resolution = 5;
t_max = 100;
t_resolution = 10;
minimum_concentration = 1e-2;

% species list
[s, species, n_species] = species_map();

reactions = [
    % photosynthesis
    1, s('C(IV)'), 0, s('photons'), 2, s('O(0)'), 1, s('C(0)'), photo_delta_G_standard, 1
    
    % denitrification
    % N(V) + 2.0C(0) -> N(-III) + 2.0C(IV): delta Go  = -3.361363e-04
    2.0, s('C(0)'), 1.0, s('N(V)'), 2.0, s('C(IV)'), 1.0, s('N(-III)'), -3.6e-4, 0
    
    % ammonia oxidation
    % O(0) + 0.25N(-III) -> water + 0.25N(V): delta Go  = 1.505088e-04
    %1, s('N(-III)'), 4, s('O(0)'), 1, s('N(V)'), 0, s('water'), -8.1e-5, 0

];
n_reactions = n_rows(reactions);

% initial conditions
function [u] = icfun(x)
    u = repmat(minimum_concentration, n_species, 1);
        
    % make water so it doesn't enter the photosynthesis ln Q
    u(s('water')) = 1.0;

    % photon density decays exponentially
    u(s('photons')) = exp(-x / photon_depth_scale);
    u(s('O(0)')) = 100 * u(s('photons'));
    
    u(s('Fe(II)')) = 100;
    
    u(s('C(IV)')) = 100;
    u(s('C(-IV)')) = 100;
    
    u(s('S(VI)')) = 100;
    
    u(s('N(V)')) = 1;
    u(s('N(-III)')) = 100;
end

% boundary conditions
function [pl, ql, pr, qr] = bcfun(xl, ul, xr, ur, t)
    pl = zeros(n_species, 1);
    ql = ones(n_species, 1);
    pr = zeros(n_species, 1);
    qr = ones(n_species, 1);
end

% computed simulation parameters
xmesh = linspace(0, x_max, x_resolution);
tspan = linspace(0, t_max, t_resolution);



function [so] = source(x, u)
    %SOURCE compute the fluxes from the concentration vector u
    so = zeros(1, length(u));
    
    for i = 1: n_reactions
        % get the reaction data from the reaction matrix
        reac1_coeff = reactions(i, 1);
        reac1_i = reactions(i, 2);
        reac1 = u(reac1_i);

        reac2_coeff = reactions(i, 3);
        reac2_i = reactions(i, 4);
        reac2 = u(reac2_i);
        
        prod1_coeff = reactions(i, 5);
        prod1_i = reactions(i, 6);
        prod1 = u(prod1_i);

        prod2_coeff = reactions(i, 7);
        prod2_i = reactions(i, 8);
        prod2 = u(prod2_i);

        delta_G_standard = reactions(i, 9);

        % compute the actual delta G   
        ln_Q = prod1_coeff * log(prod1) + prod2_coeff * log(prod2) - reac1_coeff * log(reac1) - reac2_coeff * log(reac2);
        delta_G = delta_G_standard + RT * ln_Q;

        % check if the reaction will proceed
        if reac1 > 0.0 && reac2 > 0.0 && delta_G < 0.0 && abs(delta_G) > abs(metabolic_cutoff)
            rate = rate_constant * abs(delta_G - metabolic_cutoff);
            assert(rate > 0)
            
            % if this is photosynthesis, also check for the number of photons
            if reactions(i, 10) == 1
                rate = rate * u(s('photons'));
            end

            so(reac1_i) = so(reac1_i) - reac1_coeff * rate;
            so(reac2_i) = so(reac2_i) - reac2_coeff * rate;
            so(prod1_i) = so(prod1_i) + prod1_coeff * rate;
            so(prod2_i) = so(prod2_i) + prod2_coeff * rate;
        end
        
        {'rxni', i, 'dg', delta_G, 'dgo', delta_G_standard, 'lnq', ln_Q, 'rate', rate, reac1, reac2, prod1, prod2, 'lnK', -delta_G_standard / RT};
    end
    
    % methanogenesis
    
    % check that everything went somewhere
    % but we pull O2 from nowhere using photosynthesis!
    %%assert(sum(so) == 0)
end



% -- function for computing the instantaneous differential equations --
function [c, f, so] = pdefun(x, t, u, dudx)

    
    % check that all the concentrations are positive
    if min(u) < 0 && 1 == 0
        [c, i] = min(u);
        {'t' t 'u' u c i}
        u
        species(i)
    end
    %assert(min(u) >= 0)    % enforce nonnegativity in a kludge way
    %assert(min(u) > 0)
    %u(u < minimum_concentration) = minimum_concentration;
    
    % coupling constant is 1 for each species
    c = ones(n_species, 1);

    % diffusion term is the same for all species except photons and water
    f = diffusion_constant * dudx;
    f(s('water')) = 0.0;
    f(s('photons')) = 0.0;

    % compute the source terms from the previously defined function
    so = source(x, u);
    
    % ignore any changes to water or photon concentration
    so(s('water')) = 0.0;
    so(s('photons')) = 0.0;
end

m = 1;

use_options = 0;
if use_options
    options = odeset('RelTol', 1e-3 * minimum_concentration, 'MaxStep', 1e-4);
    sol = pdepe(m, @pdefun, @icfun, @bcfun, xmesh, tspan, options);
else
    sol = pdepe(m, @pdefun, @icfun, @bcfun, xmesh, tspan);
end

end
