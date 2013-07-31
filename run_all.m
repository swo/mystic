s = species_map();
sol = run();

clf;

hold on;
for i = [s('O(0)'), s('C(0)'), s('C(IV)')]
    u = sol(:, :, i);
    u(u < 0) = 1e-7;
    surf(u)
end
xlabel('depth');
ylabel('time');
zlabel('concentration');
set(gca, 'ZScale', 'log');
hold off;