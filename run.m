function [t, y] = run()

% constants
diffusion_constant = 1.0;
precipitation_constant = 1.0 * diffusion_constant;
rate_constant = 1e-3;
max_rate = 1e12;

% simulation parameters
n_x = 6;
t_max = 1e3;

% species list
[s, ~, n_species] = species_map();

biotic_rxns = [
    % respiration
    %s('C'), s('O'), s(''), s(''), -10.0
    s('C'), s('N+'), s('N-'), s(''), -350.0
    s('C'), s('Fe+'), s('Fe-'), s(''), -250.0
    s('C'), s('S+'), s('S-'), s(''), -40.0

    % oxidations
    s('O'), s('N-'), s('N+'), s(''), -300.0
    s('O'), s('S-'), s('S+'), s(''), -300.0
    
    % fermentation
    %s('C'), s(''), s(''), s(''), -1.0
];

abiotic_rxns = [
    % iron oxidation
    s('O'), s('Fe-'), s('Fe+'), s(''), -1e4
    %s(''), s(''), s(''), s(''), -200.0
];

photo = 1e5;
sources = [
    s('O'), photo, 1
    s('C'), photo, 1
];

precipitating_species = [s('Fe+')];
%precipitating_species = [];

diffusing_species = setdiff(1: n_species, precipitating_species);

% initialize the lake
concs0 = zeros(n_x, n_species);
concs0(:, s('')) = 1.0;
concs0(:, s('C')) = 0.0;
concs0(:, s('O')) = 0.0;

concs0(:, s('N+')) = 0.0;
concs0(:, s('N-')) = 100.0;

concs0(:, s('Fe+')) = 0.0;
concs0(:, s('Fe-')) = 100.0;

concs0(:, s('S+')) = 0.0;
concs0(:, s('S-')) = 100.0;


% compute internal parameters

% grab the unchanging columns from the reaction matrix
% delta_Go is transposed so it can fit with row vectors of concentrations
bio_reac1_i = biotic_rxns(:, 1);
bio_reac2_i = biotic_rxns(:, 2);
bio_prod1_i = biotic_rxns(:, 3);
bio_prod2_i = biotic_rxns(:, 4);
bio_delta_Go = biotic_rxns(:, 5)';

abio_reac1_i = abiotic_rxns(:, 1);
abio_reac2_i = abiotic_rxns(:, 2);
abio_prod1_i = abiotic_rxns(:, 3);
abio_prod2_i = abiotic_rxns(:, 4);
abio_delta_Go = abiotic_rxns(:, 5)';


source_species = sources(:, 1);
source_rates = sources(:, 2);
source_x = sources(:, 3);
source_idx = sub2ind([n_x, n_species], source_x, source_species);


% -- flux function --
n_total = n_x * n_species;
function [fluxes] = flux(~, concs_vector)
    concs = reshape(concs_vector, [n_x, n_species]);
    fluxes = zeros(n_x, n_species);
    
    % apply the sources
    fluxes(source_idx) = fluxes(source_idx) + source_rates;

    for x = 1: n_x
        % -- reactions --
        % column column vectors of concentrations
        bio_reac1 = concs(x, bio_reac1_i);
        bio_reac2 = concs(x, bio_reac2_i);

        abio_reac1 = concs(x, abio_reac1_i);
        abio_reac2 = concs(x, abio_reac2_i);

        bio_rates = max(0.0, rate_constant * bio_reac1 .* bio_reac2 .* (-bio_delta_Go));
        abio_rates = rate_constant * abio_reac1 .* abio_reac2 .* (-abio_delta_Go);
        
        if any(bio_rates < 0.0)
            concs
            bio_rates
        end
        %assert(all(bio_rates >= -1e-6))

        % enforce the maximum rate for biotic reactions
        if sum(bio_rates) > max_rate
            %'over max rate'
            bio_rates = bio_rates * max_rate / sum(bio_rates);
        end
        
        % add rates to species
        fluxes(x, :) = fluxes(x, :) - accumarray(bio_reac1_i, bio_rates, [n_species, 1])';
        fluxes(x, :) = fluxes(x, :) - accumarray(bio_reac2_i, bio_rates, [n_species, 1])';
        fluxes(x, :) = fluxes(x, :) + accumarray(bio_prod1_i, bio_rates, [n_species, 1])';
        fluxes(x, :) = fluxes(x, :) + accumarray(bio_prod2_i, bio_rates, [n_species, 1])';
        
        fluxes(x, :) = fluxes(x, :) - accumarray(abio_reac1_i, abio_rates, [n_species, 1])';
        fluxes(x, :) = fluxes(x, :) - accumarray(abio_reac2_i, abio_rates, [n_species, 1])';
        fluxes(x, :) = fluxes(x, :) + accumarray(abio_prod1_i, abio_rates, [n_species, 1])';
        fluxes(x, :) = fluxes(x, :) + accumarray(abio_prod2_i, abio_rates, [n_species, 1])';

        % -- diffusion --
        if x > 1
            fluxes(x, diffusing_species) = fluxes(x, diffusing_species) + diffusion_constant * (concs(x - 1, diffusing_species) - concs(x, diffusing_species));
            fluxes(x, precipitating_species) = fluxes(x, precipitating_species) + precipitation_constant * concs(x - 1, precipitating_species);
        end

        if x < n_x
            fluxes(x, diffusing_species) = fluxes(x, diffusing_species) + diffusion_constant * (concs(x + 1, diffusing_species) - concs(x, diffusing_species));
        end

    end % for x
    
    % nullify fluxes on null species
    fluxes(:, s('')) = 0.0;
    
    
    fluxes = reshape(fluxes, [n_total, 1]);
end


% -- run the ode solver --
options = odeset('NonNegative', 1: n_total);
concs0_vector = reshape(concs0, [n_total, 1]);
[t, y] = ode15s(@flux, [0 t_max], concs0_vector, options);

[n_time_slices, ~] = size(y);
y = reshape(y, n_time_slices, n_x, n_species);

end
