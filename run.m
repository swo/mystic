function [sol] = run()

% constants
diffusion_constant = 0.0;
RT = 2.49; % 8.3 J K^-1 mol^-1 * 300 K = 2.49 kJ mol^-1
rate_constant = 1.0;
faraday_constant = 96.5; % 96.4 kJ volt^-1 mol^-1

% assertive parameters
photon_depth_scale = 0.5; % 1/e distance for photosynthesis (meters)
photo_delta_G_standard = -500;

% methanogenesis parameters
mg_rate_constant = 1.0;
mg_delta_G_modifier = 0.0; % accounts for the [H20]^2/[H2]^4 in Q

metabolic_cutoff = 0.0; % Canfield's cutoff for useful metabolism; -20 kJ mol^-1
assert(metabolic_cutoff <= 0.0)

% simulation parameters
x_max = 15;
x_resolution = 10;
t_max = 100;
t_resolution = 10;
minimum_concentration = 1e-9;

% species list
[s, species, n_species] = species_map();

reactions = [
    % photosynthesis
    1, s('C(IV)'), 0, s('water'), 2, s('O(0)'), 1, s('C(0)'), photo_delta_G_standard
    
    % denitrification
    % N(V) + 2.0C(0) -> N(-III) + 2.0C(IV): delta Go  = -3.361363e-04
    2.0, s('C(0)'), 1.0, s('N(V)'), 2.0, s('C(IV)'), 1.0, s('N(-III)'), -364
    
    % ammonia oxidation
    % O(0) + 0.25N(-III) -> water + 0.25N(V): delta Go  = 1.505088e-04
    1, s('N(-III)'), 4, s('O(0)'), 1, s('N(V)'), 1, s('water'), -81

];
n_reactions = n_rows(reactions);
photosynthesis_i = 1;   % index of photosynthesis reaction

% initial conditions
function [u] = icfun(x)
    u = repmat(minimum_concentration, n_species, 1);
        
    % make water so it doesn't enter the photosynthesis ln Q
    u(s('water')) = 1.0;

    % photon density decays exponentially
    %u(s('photons')) = exp(-x / photon_depth_scale);
    %u(s('O(0)')) = 10e-6 * u(s('photons'));
    u(s('O(0)')) = 10e-6;
    
    u(s('Fe(II)')) = 150e-6;
    u(s('Fe(III)')) = 10e-6;
    
    u(s('C(IV)')) = 1500e-6;
    u(s('C(-IV)')) = 10e-6;
    
    u(s('S(VI)')) = 100e-6;
    u(s('S(-II)')) = 10e-6;
    
    u(s('N(V)')) = 10e-6;
    u(s('N(-III)')) = 100e-6;
    
    % convert to log domain
    u = log10(u);
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

    % column column vectors of concentrations
    reac1 = u(reac1_i);
    reac2 = u(reac2_i);
    prod1 = u(prod1_i);
    prod2 = u(prod2_i);
    
    % compute the actual delta G
    ln_Q = prod1_coeff .* prod1 + prod2_coeff .* prod2 - reac1_coeff .* reac1 - reac2_coeff .* reac2;
    if any(isnan(ln_Q))
        u
        prod1_coeff
        prod1
        prod2_coeff
        prod2
        reac1_coeff
        reac1
        reac2_coeff
        reac2
    end
    assert(~any(isnan(ln_Q)))
    delta_G = delta_G_standard + RT * ln_Q;
    assert(~any(isnan(delta_G)))

    % ignore reactions that go backwards
    delta_G(delta_G > 0.0) = 0.0;
    if ~all(delta_G <= 0.0)
        delta_G
    end
    assert(all(delta_G <= 0.0))

    % decrease all reactions by the metabolic cutoff
    %rate = rate_constant * reac1 .^ reac1_coeff .* reac2 .^ reac2_coeff .* max(0, -delta_G + metabolic_cutoff);
    %rate = -rate_constant * delta_G;
    rate = rate_constant * (reac1_coeff .* exp10(reac1) + reac2_coeff .* exp10(reac1)) .* (-delta_G)
            
    % if this is photosynthesis, also check for the number of photons
    rate(photosynthesis_i) = rate(photosynthesis_i) * exp(-x / photon_depth_scale)
    assert(all(rate >= 0.0))
    assert(all(isfinite(rate)))

    % add rates to species, converting to log domains as they go by
    % dividing by the concentrations
    so(reac1_i) = so(reac1_i) - reac1_coeff .* rate;
    so(reac2_i) = so(reac2_i) - reac2_coeff .* rate;
    so(prod1_i) = so(prod1_i) + prod1_coeff .* rate;
    so(prod2_i) = so(prod2_i) + prod2_coeff .* rate;
    
    so = so'
end


% -- function for computing the instantaneous differential equations --
ignored_species = [s('water')];
function [c, f, so] = pdefun(x, t, u, dudx)
    if ~all(isfinite(dudx))
        u
        dudx
        t
    end
    assert(all(isfinite(dudx)))
    assert(all(isfinite(u)))
    assert(~any(isnan(u)))
    assert(~any(isnan(dudx)))
    
    % coupling constant is 1 for each species
    c = ones(n_species, 1);

    % diffusion term is the same for all species except photons and water
    f = diffusion_constant * dudx;
    f(ignored_species) = 0.0;
    assert(all(isfinite(f)))
    
    % add the second half of the diffusion terms using the source term
    so = diffusion_constant * dudx' .^ 2;
    so(ignored_species) = 0.0;

    % compute the source terms from the previously defined function
    'from before'
    so
    'exp10'
    exp10(u')
    'after source and exp added'
    so = so + source(x, u) ./ exp10(u')
    'done adding'
    assert(~any(isnan(so)))
    assert(all(isfinite(so)))
    assert(all(abs(so) < 10))
    
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
