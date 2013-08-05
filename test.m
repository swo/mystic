% todo:
% add an oxygen production term: oxygen can only be produced near the
% surface, but it requires carbon species that flow from the top

% think again about the carbon cycle

function [sol] = test()

% constants
diffusion_constant = 0.0;
RT = 1.0;
rate_constant = 1e-1;

% simulation parameters
x_max = 15;
x_resolution = 9;
t_max = 100;
t_resolution = 10;

% species list
species = {
    'A',
    'B',
    'a',
    'b'
};
n_species = length(species);
s = containers.Map(species, 1: n_species);

reactions = [
    1, s('A'), 1, s('B'), 1, s('a'), 1, s('b'), -2.0

];
n_reactions = n_rows(reactions);

% initial conditions
function [u] = icfun(x)
    u = [2, 2, 1, 1];
    
    % convert to log domain
    %u = log(u);
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
    so = zeros(1, length(u));

    % column column vectors of concentrations
    reac1 = u(reac1_i);
    reac2 = u(reac2_i);
    prod1 = u(prod1_i);
    prod2 = u(prod2_i);
    
    % compute the actual delta G
    ln_Q = prod1_coeff .* prod1 + prod2_coeff .* prod2 - reac1_coeff .* reac1 - reac2_coeff .* reac2;
    delta_G = delta_G_standard + RT * ln_Q;

    % ignore reactions that go backwards
    delta_G(delta_G > 0.0) = 0.0;

    % decrease all reactions by the metabolic cutoff
    %rate = rate_constant * reac1 .^ reac1_coeff .* reac2 .^ reac2_coeff .* max(0, -delta_G + metabolic_cutoff);
    rate = -rate_constant * delta_G;
    assert(all(rate >= 0.0))

    % add rates to species, converting to log domains as they go by
    % dividing by the concentrations
    so(reac1_i) = so(reac1_i) - reac1_coeff .* rate ./ exp(reac1);
    so(reac2_i) = so(reac2_i) - reac2_coeff .* rate ./ exp(reac2);
    so(prod1_i) = so(prod1_i) + prod1_coeff .* rate ./ exp(prod1);
    so(prod2_i) = so(prod2_i) + prod2_coeff .* rate ./ exp(prod2);
    
end


% -- function for computing the instantaneous differential equations --
function [c, f, so] = pdefun(x, t, u, dudx)
    
    % coupling constant is 1 for each species
    c = ones(n_species, 1);

    % diffusion term is the same for all species except photons and water
    f = diffusion_constant * dudx;
    
    % add the second half of the diffusion terms using the source term
    so = diffusion_constant * dudx' .^ 2;

    % compute the source terms from the previously defined function
    so = so + source(x, u);

end

% m=0 -> slab geometry, no cylindrical or spherical
m = 0;
sol = pdepe(m, @pdefun, @icfun, @bcfun, xmesh, tspan);

end