% todo:
% add an oxygen production term: oxygen can only be produced near the
% surface, but it requires carbon species that flow from the top

% think again about the carbon cycle

function [sol] = run()

% constants
diffusion_constant = 0.1;
RT = 2.49e-6; % 8.3 J K^-1 mol^-1 * 300 K = 2.49 kJ mol^-1 = 2.49e-6 kJ umol^-1
rate_constant = 1e1;
faraday_constant = 9.6485e-5; % 96.4 kJ volt^-1 mol^-1 = 96.4e-6 umol^-1

% assertive parameters
photon_depth_scale = 1.0; % 1/e distance for photosynthesis (meters)
photo_delta_G_standard = -1e-3;

% methanogenesis parameters
mg_rate_constant = 1.0;
mg_delta_G_modifier = 0.0; % accounts for the [H20]^2/[H2]^4 in Q

metabolic_cutoff = 0.0; % Canfield's cutoff for useful metabolism; -20 kJ mol^-1 = -2e-5 kJ mmol^-1
assert(metabolic_cutoff <= 0.0)

% simulation parameters
x_max = 15;
x_resolution = 5;
t_max = 1000;
t_resolution = 10;
minimum_concentration = 1e-2;

% species list
[s, species, n_species] = species_map();

reactions = [
    % photosynthesis
    1, s('C(IV)'), 0, s('photons'), 2, s('O(0)'), 1, s('C(0)'), photo_delta_G_standard
    
    % denitrification
    % N(V) + 2.0C(0) -> N(-III) + 2.0C(IV): delta Go  = -3.361363e-04
    2.0, s('C(0)'), 1.0, s('N(V)'), 2.0, s('C(IV)'), 1.0, s('N(-III)'), -3.6e-4
    
    % ammonia oxidation
    % O(0) + 0.25N(-III) -> water + 0.25N(V): delta Go  = 1.505088e-04
    1, s('N(-III)'), 4, s('O(0)'), 1, s('N(V)'), 1, s('water'), -8.1e-5

];
n_reactions = n_rows(reactions);
photosynthesis_i = 1;   % index of photosynthesis reaction

% initial conditions
function [u] = icfun(x)
    u = repmat(minimum_concentration, n_species, 1);
        
    % make water so it doesn't enter the photosynthesis ln Q
    u(s('water')) = 1.0;

    % photon density decays exponentially
    u(s('photons')) = exp(-x / photon_depth_scale);
    u(s('O(0)')) = 10 * u(s('photons'));
    
    u(s('Fe(II)')) = 150;
    
    u(s('C(IV)')) = 1500;
    u(s('C(-IV)')) = 100;
    
    u(s('S(VI)')) = 100;
    
    u(s('N(V)')) = 50;
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


% -- source term --
% grab the unchanging columns from the reaction matrix
reac1_coeff = reactions(:, 1);
reac1_i = reactions(:, 2);
reac2_coeff = reactions(:, 3);
reac2_i = reactions(:, 4);

prod1_coeff = reactions(:, 5);
prod1_i = reactions(:, 6);
prod2_coeff = reactions(:, 7);
prod2_i = reactions(:, 8);

delta_G_standard = reactions(:, 9);

function [so] = source(x, u)
    %SOURCE compute the fluxes from the concentration vector u
    so = zeros(length(u), 1);
    u
    u(s('C(0)'))
    
    if any(u < 0)
        species(u < 0)
    end
    assert(all(u > 0))

    % column column vectors of concentrations
    reac1 = u(reac1_i);
    reac2 = u(reac2_i);
    prod1 = u(prod1_i);
    prod2 = u(prod2_i);
    
    % compute the actual delta G
    ln_Q = prod1_coeff .* log(prod1) + prod2_coeff .* log(prod2) - ...
        (reac1_coeff .* log(reac1) + reac2_coeff .* log(reac2));
    delta_G = delta_G_standard + RT * ln_Q;
    
    assert(max(imag(ln_Q)) == 0)

    % ignore reactions that go backwards
    delta_G(delta_G > 0.0) = 0.0;

    % stop reactions that have negative reactants
    delta_G(reac1 < 0.0 | reac2 < 0.0) = 0.0;

    % decrease all reactions by the metabolic cutoff
    rate = rate_constant * reac1 .^ reac1_coeff .* reac2 .^ reac2_coeff .* max(0, -delta_G + metabolic_cutoff);
            
    % if this is photosynthesis, also check for the number of photons
    rate(photosynthesis_i) = rate(photosynthesis_i) * u(s('photons'));

    so(reac1_i) = so(reac1_i) - reac1_coeff .* rate;
    so(reac2_i) = so(reac2_i) - reac2_coeff .* rate;
    so(prod1_i) = so(prod1_i) + prod1_coeff .* rate;
    so(prod2_i) = so(prod2_i) + prod2_coeff .* rate;
    
    so = so.'
end


% -- function for computing the instantaneous differential equations --
ignored_species = [s('water') s('photons')];
function [c, f, so] = pdefun(x, t, u, dudx)
    t
    
    % coupling constant is 1 for each species
    c = ones(n_species, 1);

    % diffusion term is the same for all species except photons and water
    f = diffusion_constant * dudx;
    f(ignored_species) = 0.0;

    % compute the source terms from the previously defined function
    so = source(x, u);
    
    % ignore any changes to water or photon concentration
    so(ignored_species) = 0.0;
end

% m=0 -> slab geometry, no cylindrical or spherical
m = 0;

use_options = 0;
if use_options
    options = odeset('RelTol', 1e-3 * minimum_concentration, 'MaxStep', 1e-4);
    sol = pdepe(m, @pdefun, @icfun, @bcfun, xmesh, tspan, options);
else
    sol = pdepe(m, @pdefun, @icfun, @bcfun, xmesh, tspan);
end

end
