function [sol] = run_all()

% constants
diffusion_constant = 0.1;
RT = 2.49; % kJ mol^-1; T = 300 K
rate_constant = 1;

% simulation parameters
x_max = 10;
x_resolution = 10;
t_max = 0.1;
t_resolution = 20;
minimum_concentration = 1e-7;

% species list
% specify the names of all the species and use those as keys that uses a
% map 's' to link them to integers. these integers will be their position
% in the reaction matrices, etc.
species = {
    'OH-',
    'O(0)',
    'Fe(III)',
    'Fe(II)',
    'S(VI)',
    'S(-II)'
};
n_species = length(species);
s = containers.Map(species, 1: n_species);

% half-reaction matrix
% rows are: (species) + (# electrons) -> (species) with (standard electrode potential)
half_reactions = [
    % oxygen
    s('O(0)')       s('OH-')    2      0.82  % Brock

    % carbon
    
    % sulfur
    s('S(VI)')      s('S(-II)') 8     0.299   % ?
    
    % iron
    s('Fe(III)')    s('Fe(II)') 1     0.769 % wikipedia table
    
    % nitrogen
    
];

% convert the last column into Gibbs free energies
% delta G = -nFE
faraday_constant = 96.485; % kJ per volt
half_reactions(:, 4) = -faraday_constant * half_reactions(:, 3) .* half_reactions(:, 4);

% keep track of the number of half reactions
n_half_reactions = n_rows(half_reactions);

function [so] = source(u)
    %SOURCE compute the fluxes from the concentration vector u
    so = zeros(1, length(u))
    
    % loop over every pair of half-reactions
    for i = 1: n_half_reactions - 1
        rxn1_reac_i = half_reactions(i, 1);
        rxn1_reac = u(rxn1_reac_i);
        rxn1_prod_i = half_reactions(i, 2);
        rxn1_prod = u(rxn1_prod_i);
        
        for j = i + 1: n_half_reactions
            rxn2_reac_i = half_reactions(j, 1);
            rxn2_reac = u(rxn2_reac_i);
            rxn2_prod_i = half_reactions(j, 2);
            rxn2_reac = u(rxn2_prod_i);
            
            % find the standard delta G
            delta_G_st = half_reactions(i, 4) - half_reactions(j, 4);
            
            % find the stoichiometric coefficient for the second reaction
            % so that the number of electrons is conserved
            coeff = half_reactions(i, 2) / half_reactions(j, 2);
            
            % compute ln Q for the forward reaction
            % NOTE! the forward reaction means that half-reaction rxn1 is
            % going forward (is a reduction) and rxn2 is going backwards (ie is an
            % oxidation)
            ln_Q = log(rxn1_prod) - log(rxn1_reac) + coeff * (log(rxn2_reac) - log(rxn2_prod));
            
            % figure out if the reaction is going to go forward or backward
            delta_G = delta_G_st + ln_Q;

            % if reaction is going backward, then swap the product and reactant pointers
            if delta_G > 0
                % swap delta G
                delta_G = -delta_G;

                % define a swap function
                swap = @(varargin) varargin{nargin:-1:1};

                % swap the indices
                [rxn1_reac_i, rxn1_prod_i] = swap(rxn1_reac_i, rxn1_prod_i);
                [rxn2_reac_i, rxn2_prod_i] = swap(rxn2_reac_i, rxn2_prod_i);

                % swap the concentrations
                [rxn1_reac, rxn1_prod] = swap(rxn1_reac, rxn1_prod);
                [rxn2_reac, rxn2_prod] = swap(rxn2_reac, rxn2_prod);

                % now we can pretend the reaction is going forward!
            end

            % reaction is going forward
            rxn1_rate = rate_constant * rxn1_reac * (rxn2_prod ^ coeff) * delta_G;
            rxn2_rate = coeff * rxn1_rate;
            
            % first half-reaction goes forward
            so(rxn1_reac_i) = so(rxn1_reac_i) - rxn1_rate;
            so(rxn1_prod_i) = so(rxn1_prod_i) + rxn1_rate;

            % second half-reaction goes backward
            so(rxn2_reac_i) = so(rxn1_reac_i) + rxn1_rate;
            so(rxn2_prod_i) = so(rxn1_prod_i) - rxn1_rate;
end
    
% -- compute the full reaction matrix --

reactions = zeros(n_reactions, n_species);
delta_G_standard = zeros(n_reactions, 1);
k = 1;
for i = 1: n_half_reactions - 1
    % keep track of the two species involved in the first half reaction
    % (reduction)
    rxn1_reac = half_reactions(i, 1);
    rxn1_prod = half_reactions(i, 3);
    for j = i + 1: n_half_reactions
        % and in the second half-reaction (oxidation)
        rxn2_reac = half_reactions(j, 1);
        rxn2_prod = half_reactions(j, 3);
        
        % get standard Gibbs free energy of reaction, reversing the sign of
        % the oxidation half-reaction
        delta_G_standard(k) = half_reactions(i, 4) - half_reactions(j, 4);
        
        % make coefficients of reduced species 1
        reactions(k, rxn1_reac) = -1;
        reactions(k, rxn1_prod) = 1;
        
        % adjust coefficients of oxidized species based on relative number
        % of electrons in the two half-reactions
        coeff = half_reactions(i, 2) / half_reactions(j, 2);
        
        % signs in the oxidation are reversed (since now the reduced
        % species is being used up)
        reactions(k, rxn2_prod) = -coeff;
        reactions(k, rxn2_reac) = coeff;
        
        % increment the counter for the reaction number
        k = k + 1;
    end
end

% initial conditions
function [u] = icfun(x)
    u = repmat(minimum_concentration, n_species, 1);
    u(s('OH-')) = 1e-7;
    u(s('O(0)')) = exp(-x);
    u(s('Fe(II)')) = exp(-(x_max - x));
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

% -- function for computing the instantaneous differential equations --
function [c, f, so] = pdefun(x, t, u, dudx)
    
    % coupling constant is 1 for each species
    c = ones(n_species, 1);

    % diffusion term is the same for all species
    f = diffusion_constant * dudx;

    % source term is determined from all the Gibbs free energies
    % delta G = delta G standard + ln Q, so compute ln_Q
    % (ln_Q)_i = sum_j reactions_i,j * (ln u)_j
    % sum(A, 2) returns a vector of the sums of each row of A
    ln_Q = sum(bsxfun(@times, reactions, (log(u)')), 2);
    
    % compute delta G's
    delta_G = delta_G_standard + RT * ln_Q;
    
    % compute the reaction rates: multiply each row in the reaction matrix
    % by the corresponding delta G, then sum up all the reactions (i.e.,
    % sum up each row; i.e., sum by column)
    so = -rate_constant * sum(bsxfun(@times, reactions, delta_G));
    
    % ignore any changes to hydroxide concentration
    so(s('OH-')) = 0;
    
    % check that all the concentrations are positive
    if min(u) < 0
        'holy shit'
        [v, i] = min(u)
        species(i)
        u
        so
        ln_Q
        delta_G_standard
    end
    assert(min(u) >= 0)
end

m = 1;
options = odeset('AbsTol', 1e-2 * minimum_concentration);
sol = pdepe(m, @pdefun, @icfun, @bcfun, xmesh, tspan, options);

end
