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
    u = ones(n_species, 1);
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

dummy, n_species = size(delta_G_standard);

m = 1;
sol = pdepe(m, pdefun, icfun, bcfun, xmesh, tspan)

end