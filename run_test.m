% todo:
% add an oxygen production term: oxygen can only be produced near the
% surface, but it requires carbon species that flow from the top

% think again about the carbon cycle

function [so] = run_test()

% constants
diffusion_constant = 0.1;
RT = 2.49e-6; % 8.3 J K^-1 mol^-1 * 300 K = 2.49 kJ mol^-1 = 2.49e-6 kJ umol^-1
rate_constant = 1e2;
faraday_constant = 9.6485e-5; % 96.4 kJ volt^-1 mol^-1 = 96.4e-6 umol^-1

% assertive parameters
photon_depth_scale = 1.0; % 1/e distance for photosynthesis (meters)
photo_delta_G_standard = -1e-5;

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
    1, s('C(IV)'), 1, s('photons'), 2, s('O(0)'), 1, s('C(0)'), photo_delta_G_standard
    
    % denitrification
    % N(V) + 2.0C(0) -> N(-III) + 2.0C(IV): delta Go  = -3.361363e-04
    2.0, s('C(0)'), 1.0, s('N(V)'), 2.0, s('C(IV)'), 1.0, s('N(-III)'), -3.6e-4
    
    % ammonia oxidation
    % O(0) + 0.25N(-III) -> water + 0.25N(V): delta Go  = 1.505088e-04
    1, s('N(-III)'), 4, s('O(0)'), 1, s('N(V)'), 0, s('water'), -8.1e-5

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
    %u(s('photons')) = 1.0;
    u(s('O(0)')) = 10 * u(s('photons'));
    
    u(s('Fe(II)')) = 10;
    
    u(s('C(IV)')) = 10;
    u(s('C(-IV)')) = 10;
    
    u(s('S(VI)')) = 10;
    
    u(s('N(V)')) = 1;
    u(s('N(-III)')) = 10;
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
    ln_Q = prod1_coeff .* log(prod1) + prod2_coeff .* log(prod2) - ...
        (reac1_coeff .* log(reac1) + reac2_coeff .* log(reac2));
    delta_G = delta_G_standard + RT * ln_Q;

    % ignore reactions that go backwards
    delta_G(delta_G > 0.0) = 0.0;

    % stop reactions that have negative reactants
    delta_G(reac1 < 0.0 | reac2 < 0.0) = 0.0;

    % decrease all reactions by the metabolic cutoff
    rate = rate_constant * max(0, -delta_G + metabolic_cutoff);
            
    % if this is photosynthesis, also check for the number of photons
    rate(photosynthesis_i) = rate(photosynthesis_i) * u(s('photons'));

    so(reac1_i) = so(reac1_i) - reac1_coeff .* rate;
    so(reac2_i) = so(reac2_i) - reac2_coeff .* rate;
    so(prod1_i) = so(prod1_i) + prod1_coeff .* rate;
    so(prod2_i) = so(prod2_i) + prod2_coeff .* rate;
    
    so = so.';
end


species = {
    'photons',
    'water',
    'O(0)',
    'Fe(III)',
    'Fe(II)',
    'C(-IV)',
    'C(0)',
    'C(IV)'
    'S(VI)',
    'S(-II)',
    'N(V)',
    'N(-III)'
};
ignored_species = [s('photons'), s('water')];

n = 100;
k = length(species);
u = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1].';
u = randsample(n, k);

x = 0.0;
so = source(x, u);
so(ignored_species) = 0.0;

end
