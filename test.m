[x, y] = meshgrid(0: 1: 100);
kp = 1.0;
km = 1.0;
c1 = 1.0;
c2 = 1.0;

z1 = kp * x - km * y;
z3 = kp * log(x) - km * log(y);
z2 = c1 + c2 * log(x ./ y);

clf
hold on
contour(x, y, z3)
%set(gca, 'ZScale', 'log')
hold off