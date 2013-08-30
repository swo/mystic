function [sol] = run()

% constants
diffusion_constant = 1.0;
RT = 2.49; % 8.3 J K^-1 mol^-1 * 300 K = 2.49 kJ mol^-1
rate_constant = 1e-3;
max_rate = 1e2;

oxygen_constant = 10.0;
carbon_constant = 10.0;

% simulation parameters
x_max = 15;
x_resolution = 10;
t_max = 10.0;
t_resolution = 50;

% species list
[s, species, n_species] = species_map();

biotic_rxns = [
    % respiration
    s('C'), s('O'), s(''), s(''), -100.0
    s('C'), s('N+'), s('N-'), s(''), -80.0
    s('C'), s('Fe+'), s('Fe-'), s(''), -60.0
    s('C'), s('S+'), s('S-'), s(''), -40.0

    % oxidations
    s('O'), s('N-'), s('N+'), s(''), -50.0
    s('O'), s('S-'), s('S+'), s(''), -50.0
];

abiotic_rxns = [
    % iron oxidation
    s('O'), s('Fe-'), s('Fe+'), s(''), -100.0
];

n_biotic_rxns = n_rows(biotic_rxns);
n_abiotic_rxns = n_rows(abiotic_rxns);

% grab the unchanging columns from the reaction matrix
bio_reac1_i = biotic_rxns(:, 1);
bio_reac2_i = biotic_rxns(:, 2);
bio_prod1_i = biotic_rxns(:, 3);
bio_prod2_i = biotic_rxns(:, 4);
bio_delta_Go = biotic_rxns(:, 5);

abio_reac1_i = abiotic_rxns(:, 1);
abio_reac2_i = abiotic_rxns(:, 2);
abio_prod1_i = abiotic_rxns(:, 3);
abio_prod2_i = abiotic_rxns(:, 4);
abio_delta_Go = abiotic_rxns(:, 5);

% initial conditions
function [u] = icfun(x)
    u = zeros(n_species, 1);

    % null species is 1.0
    u(s('')) = 1.0;

    u(s('C')) = 10.0;
    u(s('O')) = 10.0;

    u(s('N+')) = 100.0;
    u(s('N-')) = 100.0;

    u(s('Fe+')) = 150.0;
    u(s('Fe-')) = 10.0;

    u(s('S+')) = 100.0;
    u(s('S-')) = 10.0;
end

% boundary conditions
oc_i = [s('O'), s('C')];
function [pl, ql, pr, qr] = bcfun(xl, ul, xr, ur, t)    
    % mostly require that there is no flux
    pl = zeros(n_species, 1);
    ql = ones(n_species, 1);
    pr = zeros(n_species, 1);
    qr = ones(n_species, 1);

    % except for oxygen and carbon, which are kept at a constant level
    pl(s('O')) = ul(s('O')) - oxygen_constant;
    pl(s('C')) = ul(s('C')) - carbon_constant;
    ql(oc_i) = 0.0;
    qr(oc_i) = 0.0;
    
    % and zero at the bottom of the lake
    pr(s('O')) = ur(s('O'));
    pr(s('C')) = ur(s('C'));
end

% computed simulation parameters
xmesh = linspace(0, x_max, x_resolution);
tspan = linspace(0, t_max, t_resolution);


% -- source term --
function [so] = source(x, u)
    %SOURCE compute the fluxes from the concentration vector u
    so = zeros(length(u), 1);

    % column column vectors of concentrations
    bio_reac1 = max(0.0, u(bio_reac1_i));
    bio_reac2 = max(0.0, u(bio_reac2_i));
    abio_reac1 = max(0.0, u(abio_reac1_i));
    abio_reac2 = max(0.0, u(abio_reac2_i));

    bio_rates = rate_constant * bio_reac1 .* bio_reac2 .* (-bio_delta_Go);
    abio_rates = rate_constant * abio_reac1 .* abio_reac2 .* (-abio_delta_Go);

    % enforce the maximum rate for biotic reactions
    if sum(bio_rates) > max_rate
        bio_rates = bio_rates * max_rate / sum(bio_rates);
    end
    
    assert(all(abio_rates >= 0.0))
    assert(all(bio_rates >= 0.0))
    
    % add rates to species
    so(bio_reac1_i) = so(bio_reac1_i) - bio_rates;
    so(bio_reac2_i) = so(bio_reac2_i) - bio_rates;
    so(bio_prod1_i) = so(bio_prod1_i) + bio_rates;
    so(bio_prod2_i) = so(bio_prod2_i) + bio_rates;

    so(abio_reac1_i) = so(abio_reac1_i) - abio_rates;
    so(abio_reac2_i) = so(abio_reac2_i) - abio_rates;
    so(abio_prod1_i) = so(abio_prod1_i) + abio_rates;
    so(abio_prod2_i) = so(abio_prod2_i) + abio_rates;
    
    so = so';
end


% -- function for computing the instantaneous differential equations --
function [c, f, so] = pdefun(x, t, u, dudx)
    %assert(all(u > 0.0))
    
    % coupling constant is 1 for each species
    c = ones(n_species, 1);

    % diffusion term is the same for all species except photons and water
    f = diffusion_constant * dudx;

    % compute the source terms from the previously defined function
    so = source(x, u);

    % keep the null species constant
    f(s('')) = 0.0;
    so(s('')) = 0.0;
end

% m=0 -> slab geometry, no cylindrical or spherical
m = 0;

sol = pdepe(m, @pdefun, @icfun, @bcfun, xmesh, tspan);

end
