G = -10.0;
T = 5.0;

P = 2.0;
R = linspace(0.0, 0.5, 1000);
y = max(0.0, -G + T * (log(R) - log(P)));
%approx = max(0.0, R - exp(G / T) * P);

approx_bad = max(0.0, T * (R * exp(-G / (2 * T)) - P * exp(G / (2 * T))));
approx = max(0.0, T * (R ./ P * exp(-G / T) - P ./ R * exp(G / T)));

plot(R, y, 'o', R, approx, R, approx_bad)

%assert(False)

clf;
P = linspace(1.0, 10.0, 1000);
R = 1.0;
y = max(0.0, -G + T * (log(R) - log(P)));
approx = max(0.0, T * (R ./ P * exp(-G / T) - P ./ R * exp(G / T)));
approx_bad = max(0.0, T * (R * exp(-G / (2 * T)) - P * exp(G / (2 * T))));


plot(P, y, 'o', P, approx, P, approx_bad)


P = linspace(0.2, 5, 100);
R = linspace(0.1, 2, 100);

[pg, rg] = meshgrid(P, R);

true = -G + T * (log(rg) - log(pg));
approx = T * (rg ./ pg * exp(-G / T) - pg ./ rg * exp(G / T));
worse = T * (rg * exp(-G / (2*T)) - pg * exp(G / (2*T)));
bad = T * (rg - exp(G / T) * pg);
best = T * (1.0 - pg ./ rg * exp(G / T));

clf;
hold on;
zlim([-10, 10])
surf(P, R, max(0.0, true))
xlabel('P')
surf(P, R, max(0.0, approx))
%surf(P, R, approx)
hold off;

assert(False)
Q = linspace(1e-1, 1e2, 1000);
y = -G - T * log(Q);
yy = T * (1.0 - Q * exp(G / T));

plot(Q, y, Q, max(0.0, yy))