sol = run();

hold on;
for i = 1:4
    surf(sol(:, :, i))
end
hold off;