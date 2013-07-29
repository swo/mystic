sol = run();

clf;

hold on;
for i = [2 6 7]
    surf(sol(:, :, i))
end
%set(gca, 'ZScale', 'log');
hold off;