sol = run();

clf;

hold on;
for i = [2 3 4]
    surf(sol(:, :, i))
end
set(gca, 'ZScale', 'log');
hold off;