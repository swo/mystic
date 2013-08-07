s = species_map();
sol = run();

clf;

l = 0;

hold on;
for i = cellfun(@(x) s(x), {'O(0)', 'N(-III)', 'N(V)'})
    u = sol(:, :, i);
    
    surf(u)
end
xlabel('depth');
ylabel('time');
zlabel('concentration');

if l
    set(gca, 'ZScale', 'log');
end

hold off;