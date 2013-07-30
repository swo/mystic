s = species_map();
sol = run();

clf;

hold on;
for i = [s('O(0)')]
    surf(sol(:, :, i))
end
xlabel('depth');
ylabel('time');
zlabel('concentration');
%set(gca, 'ZScale', 'log');
hold off;