function [ c, f, s ] = pdefun( x, t, u, dudx )
%PDEFUN dynamic variables to next timestep

% constants
diffusion_constant = 1;
RT = 1;


[dummy, n] = size(u);

% coupling constant is 1 for each species
c = ones(1, n);

% diffusion term is the same for all species
f = diffusion_constant * dudx;

% source term is determined from all the Gibbs free energies
s = zeros(1, n);
for substrate = 1: n - 1
    for product = substrate + 1: n
        % rate of substrate -> product
        rate = delta_G_standard(substrate, product) + RT * log(u(product) / u(substrate))
        s(substrate) = s(substrate) + rate
        s(product) = s(product) - rate
    end
end
    

end

