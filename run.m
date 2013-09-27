%function [t, c, m, flux_out, bio_rates_out, abio_rates_out] = run()
function [t, c, m] = run()

% constants
diffusion_constant_per_compartment2 = 0.1;
kappa = 1.0;
mu_decay = 1e2;
efficiency = 1.0;
rate_constant = 1.0;

microbes0_value = 1.0;

fixed_oxygen_level = 100.0;
fixed_oxygen_diffusion = 1e7;

% simulation parameters
n_x = 5;
t_max = 1;
diffusion_constant = diffusion_constant_per_compartment2 * n_x ^ 2;

% species list
[s, ~, n_species] = species_map();

precipitation_constant_input = [
    s('N-'), 0.0
    s('S-'), 0.0
    s('Fe+'), 0.1
];

biotic_rxns = [
    % respiration
    s('C'), s('N+'), s('N-'), s(''), -450.0 % Canfield Table 3.7, p.90, x5/4
    s('C'), s('Fe+'), s('Fe-'), s(''), -300.0   % Ehrlich 13.6.6, p. 314
    s('C'), s('S+'), s('S-'), s(''), -80.0  % Canfield Table 3.7, x2 for 1/2

    % oxidations
    s('O'), s('N-'), s('N+'), s(''), -300.0    % Ehrlich 13.2.3, p. 235
    s('O'), s('S-'), s('S+'), s(''), -190.0 % cooky website
    
    % iron oxidation on nitrate
    s('N+'), s('Fe-'), s('Fe+'), s('N-'), -300.0    % ?? Canfield p.284 describes reaction
];

abiotic_rxns = [
    % iron oxidation
    s('O'), s('Fe-'), s('Fe+'), s(''), 1e3
];

photo = 0.0;
sources = [
    s('O'), photo, 1
    s('C'), photo, 1
];

% initialize the bacterial species
[n_microbes, ~] = size(biotic_rxns);
microbes0 = repmat(microbes0_value, n_x, n_microbes);

% initialize the lake
concs0 = zeros(n_x, n_species);
concs0(:, s('')) = 1.0;
concs0(:, s('C')) = 0.0;
concs0(:, s('O')) = 0.0;

concs0(:, s('N+')) = 100.0;
concs0(:, s('N-')) = 0.0;

concs0(:, s('Fe+')) = 30.0;
concs0(:, s('Fe-')) = 30.0;

concs0(:, s('S+')) = 60.0;
concs0(:, s('S-')) = 60.0;


% compute internal parameters

% make the precipitation constants
precipitation_constants = zeros([1 n_species]);
for i = 1: length(precipitation_constant_input)
    idx = precipitation_constant_input(i, 1);
    val = precipitation_constant_input(i, 2);

    precipitation_constants(idx) = val;
end

% make new precipitation values
D = diffusion_constant;
D_plus = (1.0 + precipitation_constants) * D;
D_minus = (1.0 - precipitation_constants) * D;

% grab the unchanging columns from the reaction matrix

% delta_Go is transposed so it can fit with row vectors of concentrations
bio_reac1_i = biotic_rxns(:, 1);
bio_reac2_i = biotic_rxns(:, 2);
bio_prod1_i = biotic_rxns(:, 3);
bio_prod2_i = biotic_rxns(:, 4);
bio_delta_Go = biotic_rxns(:, 5)';
n_bio_rxns = length(bio_delta_Go);
assert(n_bio_rxns == n_microbes)

abio_reac1_i = abiotic_rxns(:, 1);
abio_reac2_i = abiotic_rxns(:, 2);
abio_prod1_i = abiotic_rxns(:, 3);
abio_prod2_i = abiotic_rxns(:, 4);
abio_rate_constants = abiotic_rxns(:, 5)';
n_abio_rxns = length(abio_rate_constants);

source_species = sources(:, 1);
source_rates = sources(:, 2);
source_x = sources(:, 3);
source_idx = sub2ind([n_x, n_species], source_x, source_species);


% -- flux function --
n_conc_total = n_x * n_species;
n_microbes_total = n_x * n_microbes;
n_total = n_conc_total + n_microbes_total;

function [bio_rates, abio_rates, micro_rates] = rates(concs_row, micro_row)
    bio_reac1 = concs_row(bio_reac1_i);
    bio_reac2 = concs_row(bio_reac2_i);

    abio_reac1 = concs_row(abio_reac1_i);
    abio_reac2 = concs_row(abio_reac2_i);

    bio_rates = kappa * bio_reac1 .* bio_reac2 .* micro_row ./ (-bio_delta_Go);
    %assert(all(bio_rates > -1e-6))

    micro_rates = efficiency * (-bio_delta_Go) .* bio_rates - mu_decay * micro_row;

    abio_rates = rate_constant * abio_reac1 .* abio_reac2 .* abio_rate_constants;
end

% twiddle because this is time-independent
function [out_flux] = flux(~, concs_vector)
    concs = concs_vector(1: n_conc_total);
    concs = reshape(concs, [n_x, n_species]);

    micros = concs_vector(n_conc_total + 1: end);
    micros = reshape(micros, [n_x, n_microbes]);

    conc_fluxes = zeros(n_x, n_species);
    micro_fluxes = zeros(n_x, n_microbes);
    
    % apply the sources
    conc_fluxes(source_idx) = conc_fluxes(source_idx) + source_rates;
    
    % apply the fixed oxygen term
    conc_fluxes(1, [s('O') s('C')]) = conc_fluxes(1, [s('O') s('C')]) + fixed_oxygen_diffusion * (fixed_oxygen_level - concs(1, s('O')));
    
    for x = 1: n_x
        [bio_rates, abio_rates, micro_rates] = rates(concs(x, :), micros(x, :));
        %if any(bio_rates < 0.0)
        %    concs
        %    bio_rates
        %end
        
        % add rates to species
        conc_fluxes(x, :) = conc_fluxes(x, :) - accumarray(bio_reac1_i, bio_rates, [n_species, 1])';
        conc_fluxes(x, :) = conc_fluxes(x, :) - accumarray(bio_reac2_i, bio_rates, [n_species, 1])';
        conc_fluxes(x, :) = conc_fluxes(x, :) + accumarray(bio_prod1_i, bio_rates, [n_species, 1])';
        conc_fluxes(x, :) = conc_fluxes(x, :) + accumarray(bio_prod2_i, bio_rates, [n_species, 1])';
        
        conc_fluxes(x, :) = conc_fluxes(x, :) - accumarray(abio_reac1_i, abio_rates, [n_species, 1])';
        conc_fluxes(x, :) = conc_fluxes(x, :) - accumarray(abio_reac2_i, abio_rates, [n_species, 1])';
        conc_fluxes(x, :) = conc_fluxes(x, :) + accumarray(abio_prod1_i, abio_rates, [n_species, 1])';
        conc_fluxes(x, :) = conc_fluxes(x, :) + accumarray(abio_prod2_i, abio_rates, [n_species, 1])';

        micro_fluxes(x, :) = micro_rates;
        
        % -- diffusion --        
        if x > 1
            conc_fluxes(x, :) = conc_fluxes(x, :) + D_plus .* concs(x - 1, :) - D_minus .* concs(x, :);
        end

        if x < n_x
            conc_fluxes(x, :) = conc_fluxes(x, :) - D_plus .* concs(x, :) + D_minus .* concs(x + 1, :);
        end

    end % for x
    
    % nullify fluxes on null species
    conc_fluxes(:, s('')) = 0.0;
    
    conc_fluxes = reshape(conc_fluxes, [n_conc_total, 1]);
    micro_fluxes = reshape(micro_fluxes, [n_microbes_total, 1]);

    out_flux = [conc_fluxes; micro_fluxes];
end


% -- run the ode solver --
options = odeset('NonNegative', 1: n_total);
concs0_vector = reshape(concs0, [n_conc_total, 1]);
microbes0_vector = reshape(microbes0, [n_microbes_total, 1]);
vector0 = [concs0_vector; microbes0_vector];
[t, y] = ode15s(@flux, [0 t_max], vector0, options);

[n_time_slices, ~] = size(y);
c = y(:, 1: n_conc_total);
m = y(:, n_conc_total + 1: end);

c = reshape(c, n_time_slices, n_x, n_species);
m = reshape(m, n_time_slices, n_x, n_microbes);

%c_vec = reshape(c(end, :, :), [n_total, 1]);
%flux_out = reshape(flux(0, c_vec), [n_x, n_species]);

%bio_rates_out = zeros(n_x, n_bio_rxns);
%abio_rates_out = zeros(n_x, n_abio_rxns);
%for x = 1: n_x
%    [bio_rates_out(x, :), abio_rates_out(x, :)] = rates(squeeze(y(end, x, :))');
%end

end
