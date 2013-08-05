% todo:
% add an oxygen production term: oxygen can only be produced near the
% surface, but it requires carbon species that flow from the top

% think again about the carbon cycle

function [sol] = run(use_log)

% simulation parameters
x_max = 5;
x_resolution = 100;
t_max = 10;
t_resolution = 100;

xmesh = linspace(-x_max, x_max, x_resolution);
tspan = linspace(0, t_max, t_resolution);

% initial conditions
function [u] = icfun(x)
    u = 1.0 / (1.0 + x ^ 2);
    if use_log
        u = log(u);
    end
end

% boundary conditions
function [pl, ql, pr, qr] = bcfun(xl, ul, xr, ur, t)
    pl = [0];
    ql = [1];
    pr = [0];
    qr = [1];
end

function [c, f, so] = pdefun(x, t, u, dudx)
    c = [1];
    if use_log
        f = dudx;
        so = dudx .^ 2;
    else
        f = dudx;
        so = [0];
    end
end

m = 0;
sol = pdepe(m, @pdefun, @icfun, @bcfun, xmesh, tspan);

end
