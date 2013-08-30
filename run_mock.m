function [sol] = run()

% constants
diffusion_constant = 1e-2;
RT = 1.0; % 8.3 J K^-1 mol^-1 * 300 K = 2.49 kJ mol^-1
rate_constant = 1e-2;

% assertive parameters
photo_rate_constant = 1e-2;
photon_depth_scale = 0.5; % 1/e distance for photosynthesis (meters)

metabolic_cutoff = 0.0; % Canfield's cutoff for useful metabolism; -20 kJ mol^-1
assert(metabolic_cutoff <= 0.0)

% simulation parameters
x_max = 15;
x_resolution = 50;
t_max = 1.2e-4;
t_resolution = 10;
minimum_concentration = 1e-9;

% species list
species = {'', 'C0', 'CI', 'CIV', 'C-IV', 'H', 'O', 'N+', 'N-'};
n_species = length(species);
s = containers.Map(species, 1: n_species);

reactions = [
    % O reduction (aerobic respiration)
    % O + C(x) -> C(x+1)
    s('O'), s('C0'), s('CI'), s(''), -500
    s('O'), s('CI'), s('CIV'), s(''), -200

    % fermentation
    % C(x) -> C(x+1) + H
    s('C0'), s(''), s('CI'), s('H'), -100
    s('CI'), s(''), s('CIV'), s('H'), -50

    % N reduction (denitrification)
    % N+ + C(x) -> N- + C(x+1)
    s('N+'), s('C0'), s('N-'), s('CI'), -300
    s('N+'), s('CI'), s('N-'), s('CIV'), -150

    % N oxidation (ammox)
    % N- + O -> N+
    s('N-'), s('O'), s('N+'), s(''), -80
    
    % methanogenesis
    % CIV + H -> C-IV
    s('CIV'), s('H'), s('C-IV'), s(''), -8
];


% initial conditions
function [u] = icfun(x)
    u = repmat(minimum_concentration, n_species, 1);
        
    % null is at standard
    u(s('')) = 1.0;
    u(s('N-')) = 10e-6;
    
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


% -- source term --
% grab the unchanging columns from the reaction matrix
reac1_i = reactions(:, 1);
reac2_i = reactions(:, 2);

prod1_i = reactions(:, 3);
prod2_i = reactions(:, 4);

delta_G_standard = reactions(:, 5);

function [so] = source(x, u)
    %SOURCE compute the fluxes from the concentration vector u

    so = zeros(length(u), 1);

    % column column vectors of concentrations
    reac1 = u(reac1_i);
    reac2 = u(reac2_i);
    prod1 = u(prod1_i);
    prod2 = u(prod2_i);
    
    % compute the actual delta G
    ln_Q = prod1 + prod2 - reac1 - reac2;
    assert(~any(isnan(ln_Q)))

    delta_G = delta_G_standard + RT * ln_Q;
    assert(~any(isnan(delta_G)))

    % ignore reactions that go backwards
    delta_G(delta_G > 0.0) = 0.0;
    assert(all(delta_G <= 0.0))

    %rate = -rate_constant * delta_G;
    rate = rate_constant * min(exp10(reac1), exp10(reac1)) .* (-delta_G);
    %rate(photosynthesis_i) = rate(photosynthesis_i) * exp(-x / photon_depth_scale); 
    
    assert(all(rate >= 0.0))
    assert(all(isfinite(rate)))
    
    % if rate exceeds the total, then rescale the whole thing
    %if sum(rate) > total_rate
        %rate = rate / sum(rate) * total_rate;
    %end

    % now take the rates to some proportion
    % take rates so they sum to one
    %rate = total_rate * delta_G / sum(delta_G);

    % add rates to species, converting to log domains as they go by
    % dividing by the concentrations
    so(reac1_i) = so(reac1_i) - rate;
    so(reac2_i) = so(reac2_i) - rate;
    so(prod1_i) = so(prod1_i) + rate;
    so(prod2_i) = so(prod2_i) + rate;

    % photosynthesis
    so([s('O') s('C0')]) = so([s('O') s('C0')]) + photo_rate_constant * exp(-x / photon_depth_scale);
    
    so = so';
end


% -- function for computing the instantaneous differential equations --
ignored_species = [s('') s('C-IV')];
function [c, f, so] = pdefun(x, t, u, dudx)
    assert(all(isfinite(dudx)))
    assert(all(isfinite(u)))
    assert(~any(isnan(u)))
    assert(~any(isnan(dudx)))
    
    eu = exp10(u);
    
    % coupling constant is 1 for each species
    c = log(10) * eu;

    % diffusion term is the same for all species except photons and water
    f = diffusion_constant * log(10) * eu .* dudx;
    f(ignored_species) = 0.0;
    assert(all(isfinite(f)))

    % compute the source terms from the previously defined function
    %'after source and exp added'
    so = source(x, u);
    %'done adding'
    t
    assert(~any(isnan(so)))
    assert(all(isfinite(so)))
    %assert(all(abs(so) < 10))
    
    % ignore any changes to water or photon concentration
    so(ignored_species) = 0.0;
end

% m=0 -> slab geometry, no cylindrical or spherical
m = 0;
xmesh = linspace(0, x_max, x_resolution);
tspan = linspace(0, t_max, t_resolution);
sol = pdepe(m, @pdefun, @icfun, @bcfun, xmesh, tspan);

end
