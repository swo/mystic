function [s, species, n_species] = species_map()

% specify the names of all the species and use those as keys that uses a
% map 's' to link them to integers. these integers will be their position
% in the reaction matrices, etc.
species = {
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
n_species = length(species);
s = containers.Map(species, 1: n_species);

end

