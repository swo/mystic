% todo:
% add an oxygen production term: oxygen can only be produced near the
% surface, but it requires carbon species that flow from the top

% think again about the carbon cycle

function [sol] = run_all()

% constants
diffusion_constant = 0.01;
RT = 2.49e-3; % 8.3 J K^-1 mol^-1 * 300 K = 2.49 kJ mol^-1 = 2.49e-3 kJ mmol^-1
rate_constant = 1.0;

% assertive parameters
photo_depth_scale = 0.1; % 1/e distance for photosynthesis (meters)
photo_rate_constant = 1.0; % convert CO2 concentration and photon density to rate
photo_delta_G_st = -100;

metabolic_cutoff = 0.0; % Canfield's cutoff for useful metabolism

% simulation parameters
x_max = 1;
x_resolution = 10;
t_max = 1;
t_resolution = 10;
minimum_concentration = 1e-4;

% species list
% specify the names of all the species and use those as keys that uses a
% map 's' to link them to integers. these integers will be their position
% in the reaction matrices, etc.
species = {
    'OH-',
    'O(0)',
    'Fe(III)',
    'Fe(II)',
    'C(-IV)',
    'C(0)',
    'C(IV)'
    %%'S(VI)',
    %%'S(-II)'
};
n_species = length(species);
s = containers.Map(species, 1: n_species);

% half-reaction matrix
% rows are: (species) + (# electrons) -> (species) with (standard electrode potential)
half_reactions = [
    % oxygen
    s('O(0)')       s('OH-')    2       0.820  % reduction of oxygen; Brock

    % carbon
    s('C(IV)')      s('C(0)')   4      -0.071 % opposite of fermentation
    s('C(IV)')      s('C(-IV)') 8       0.170 % methanogenesis
    
    % sulfur
    %%s('S(VI)')      s('S(-II)') 8     0.299   % ?
    
    % iron
    s('Fe(III)')    s('Fe(II)') 1     0.769 % wikipedia table
    
    % nitrogen
    
];

% initial conditions
function [u] = icfun(x)
    u = repmat(minimum_concentration, n_species, 1);
    
    u(s('OH-')) = 1e-4;
    %u(s('O(0)')) = exp(-x);
    %u(s('O(0)')) = 0.1;
    u(s('Fe(II)')) = 1;
    u(s('C(IV)')) = 1;
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

% convert the last column into Gibbs free energies
% delta G = -nFE
faraday_constant = 96.485e-1; % kJ volt^-1 mmol^-1
half_reactions(:, 4) = -faraday_constant * half_reactions(:, 3) .* half_reactions(:, 4);

% keep track of the number of half reactions
n_half_reactions = n_rows(half_reactions);

function [so] = source(x, u)
    %SOURCE compute the fluxes from the concentration vector u
    so = zeros(1, length(u));
    
    % special case: photosynthesis consumes oxidized carbon (CO2) and
    % produces C(0) (glucose, carbohydrates) and oxidized oxygen (O2)
    photons = exp(-x / photo_depth_scale);
    ln_Q = log(u(s('C(0)'))) + log(u(s('O(0)'))) - log(u(s('C(IV)'))) - log(photons);
    delta_G = photo_delta_G_st + ln_Q;
    
    if abs(delta_G) > abs(metabolic_cutoff)
        photo_rate = -photo_rate_constant * u(s('C(IV)')) * photons * delta_G;
        
        %%photo_rate = photo_rate_constant * exp(-x / photo_depth_scale) * abs(minimum_concentration - u(s('C(IV)')));
        so(s('C(IV)')) = so(s('C(IV)')) - photo_rate;
        so(s('C(0)')) = so(s('C(0)')) + photo_rate;
        so(s('O(0)')) = so(s('O(0)')) + photo_rate;
    end
    
    % loop over every pair of half-reactions
    for i = 1: n_half_reactions - 1
        for j = i + 1: n_half_reactions
            % assign all these pointers inside the inner loop so they can
            % be reassigned for reversed reactions
            rxn1_reac_i = half_reactions(i, 1);
            rxn1_reac = u(rxn1_reac_i);
            rxn1_prod_i = half_reactions(i, 2);
            rxn1_prod = u(rxn1_prod_i);
            
            rxn2_reac_i = half_reactions(j, 1);
            rxn2_reac = u(rxn2_reac_i);
            rxn2_prod_i = half_reactions(j, 2);
            rxn2_prod = u(rxn2_prod_i);
            
            % find the standard delta G
            delta_G_st = half_reactions(i, 4) - half_reactions(j, 4);
            
            % find the stoichiometric coefficient for the second reaction
            % so that the number of electrons is conserved
            coeff = half_reactions(i, 3) / half_reactions(j, 3);
            
            % compute ln Q for the forward reaction
            % NOTE! the forward reaction means that half-reaction rxn1 is
            % going forward (is a reduction) and rxn2 is going backwards (ie is an
            % oxidation)
            ln_Q = log(rxn1_prod) - log(rxn1_reac) + coeff * (log(rxn2_reac) - log(rxn2_prod));
            assert(imag(ln_Q) == 0)
            
            % figure out if the reaction is going to go forward or backward
            delta_G = delta_G_st + RT * ln_Q;

            % if reaction is going backward, then swap the product and reactant pointers
            if delta_G > 0
                % swap delta G
                delta_G = -delta_G;

                % swap the indices
                [rxn1_reac_i, rxn1_prod_i] = swap(rxn1_reac_i, rxn1_prod_i);
                [rxn2_reac_i, rxn2_prod_i] = swap(rxn2_reac_i, rxn2_prod_i);

                % swap the concentrations
                [rxn1_reac, rxn1_prod] = swap(rxn1_reac, rxn1_prod);
                [rxn2_reac, rxn2_prod] = swap(rxn2_reac, rxn2_prod);
            end
            % now we can pretend the reaction is going forward!
            assert(delta_G < 0)

            % reaction is going forward
            % rate is proportional to the difference between delta G and
            % the metabolic cutoff
            if abs(delta_G) > abs(metabolic_cutoff)
                rxn1_rate = rate_constant * rxn1_reac * (rxn2_prod ^ coeff) * abs(delta_G - metabolic_cutoff);
            else
                continue
            end
            
            assert(rxn1_rate > 0)
            rxn2_rate = coeff * rxn1_rate;
            
            % first half-reaction goes forward
            so(rxn1_reac_i) = so(rxn1_reac_i) - rxn1_rate;
            so(rxn1_prod_i) = so(rxn1_prod_i) + rxn1_rate;

            % second half-reaction goes backward
            so(rxn2_reac_i) = so(rxn2_reac_i) + rxn2_rate;
            so(rxn2_prod_i) = so(rxn2_prod_i) - rxn2_rate;
        end
    end
    
    % check that everything went somewhere
    % but we pull O2 from nowhere using photosynthesis!
    %%assert(sum(so) == 0)
end



% -- function for computing the instantaneous differential equations --
function [c, f, so] = pdefun(x, t, u, dudx)
    % check that all the concentrations are positive
    if min(u) < 0
        [c i] = min(u);
        {'t' t 'u' u c i}
        u
        species(i)
    end
    assert(min(u) >= 0)
    
    % coupling constant is 1 for each species
    c = ones(n_species, 1);

    % diffusion term is the same for all species
    f = diffusion_constant * dudx;

    % compute the source terms from the previously defined function
    so = source(x, u);
    
    % ignore any changes to hydroxide concentration
    so(s('OH-')) = 0;
end

m = 1;
options = odeset('RelTol', 1e-3 * minimum_concentration, 'MaxStep', 1e-3, 'NonNegative', 7);
sol = pdepe(m, @pdefun, @icfun, @bcfun, xmesh, tspan, options);

end
