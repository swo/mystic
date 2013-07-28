sol = run();

clf;

hold on;
for i = 1:4
    surf(sol(:, :, i))
end
hold off;