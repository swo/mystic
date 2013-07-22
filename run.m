function [] = run_all()

% constants
diffusion_constant = 1;
RT = 1;
delta_G_standard = [ 0 ];
rate_constant = 1;

% simulation parameters
x_max = 10;
x_resolution = 100;
t_max = 10;
t_resolution = 100;

% initial conditions
function [u] = icfun(x)
    for i = 1: size(x)
        if 4 < x(i) < 6
            u(i) = 1;
        else
            u(i) = 0;
        end
    end
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

[dummy, n_species] = size(delta_G_standard);

function [c, f, s] = pdefun(x, t, u, dudx)
    % coupling constant is 1 for each species
    c = ones(n_species, 1);

    % diffusion term is the same for all species
    f = diffusion_constant * dudx;

    % source term is determined from all the Gibbs free energies
    s = zeros(1, n_species);
    for substrate = 1: n_species - 1
        for product = substrate + 1: n_species
            % rate of substrate -> product
            rate = rate_constant * (delta_G_standard(substrate, product) + RT * log(u(product) / u(substrate)));
            s(substrate) = s(substrate) + rate;
            s(product) = s(product) - rate;
        end
    end
   
end

m = 1;
sol = pdepe(m, @pdefun, @icfun, @bcfun, xmesh, tspan);

end