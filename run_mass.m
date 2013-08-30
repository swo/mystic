% rates go in the proportion according to their dG

function [sol] = run()

% constants
diffusion_constant = 1e3;
RT = 2.49e-6; % 8.3 J K^-1 mol^-1 * 300 K = 2.49 kJ mol^-1
rate_constant = 1.0;
total_rate = 1e-3;

% assertive parameters
photo_depth_scale = 0.5; % 1/e distance for photosynthesis (meters)
photo_rate = 1.0;
photo_delta_G_standard = -100;

% methanogenesis parameters
mg_rate_constant = 1.0;
mg_delta_G_modifier = 0.0; % accounts for the [H20]^2/[H2]^4 in Q

metabolic_cutoff = 0.0; % Canfield's cutoff for useful metabolism; -20 kJ mol^-1
assert(metabolic_cutoff <= 0.0)

% simulation parameters
x_max = 15;
x_resolution = 20;
t_max = 1000.0;
t_resolution = 10;
minimum_concentration = 1e-3;

% species list
[s, species, n_species] = species_map();

reactions = [
    % photosynthesis
    
    % denitrification
    % N(V) + 2.0C(0) -> N(-III) + 2.0C(IV): delta Go  = -3.361363e-04
    2.0, s('C(0)'), 1.0, s('N(V)'), 2.0, s('C(IV)'), 1.0, s('N(-III)'), -364
    
    % ammonia oxidation
    % O(0) + 0.25N(-III) -> water + 0.25N(V): delta Go  = 1.505088e-04
    1, s('N(-III)'), 4, s('O(0)'), 1, s('N(V)'), 1, s('water'), -81
    
    % methanogenesis
    % C(IV) -> C(-IV)
    1, s('C(IV)'), 0, s('water'), 1, s('C(-IV)'), 0, s('water'), -8.3

];
n_reactions = n_rows(reactions);
photosynthesis_i = 1;   % index of photosynthesis reaction

% initial conditions
function [u] = icfun(x)
    u = repmat(minimum_concentration, n_species, 1);
        
    % make water so it doesn't enter the photosynthesis ln Q
    u(s('water')) = 1e6;

    % photon density decays exponentially
    u(s('O(0)')) = 1;
    
    u(s('Fe(II)')) = 150;
    u(s('Fe(III)')) = 10;
    
    u(s('C(IV)')) = 1500;
    u(s('C(-IV)')) = 10;
    
    u(s('S(VI)')) = 100;
    u(s('S(-II)')) = 10;
    
    u(s('N(V)')) = 100;
    u(s('N(-III)')) = 1;
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

delta_Go = reactions(:, 9);

function [so] = source(x, u)
    %SOURCE compute the fluxes from the concentration vector u
    so = zeros(length(u), 1);
    
    % add photosynthesis
    so(s('O(0)')) = photo_rate * exp(-x / photo_depth_scale);
    so(s('C(0)')) = photo_rate * exp(-x / photo_depth_scale);

    % column column vectors of concentrations
    reac1 = u(reac1_i);
    reac2 = u(reac2_i);
    prod1 = u(prod1_i);
    prod2 = u(prod2_i);
    
    % compute rates
    rate = rate_constant * RT * max(0, reac1 .* reac2 - exp(delta_Go / RT) .* prod1 .* prod2);

    % decrease all reactions by the metabolic cutoff

    assert(all(rate >= 0.0))
    assert(all(isfinite(rate)))
    
    % if rate exceeds the total, then rescale the whole thing

    % add rates to species
    so(reac1_i) = so(reac1_i) - reac1_coeff .* rate;
    so(reac2_i) = so(reac2_i) - reac2_coeff .* rate;
    so(prod1_i) = so(prod1_i) + prod1_coeff .* rate;
    so(prod2_i) = so(prod2_i) + prod2_coeff .* rate;
    
    so = so';
end


% -- function for computing the instantaneous differential equations --
ignored_species = [s('water')];
function [c, f, so] = pdefun(x, t, u, dudx)
    
    % coupling constant is 1 for each species
    c = ones(n_species);

    % diffusion term is the same for all species except photons and water
    f = diffusion_constant * dudx;
    f(ignored_species) = 0.0;
    assert(all(isfinite(f)))

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
