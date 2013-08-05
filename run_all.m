s = species_map();
sol = run();

clf;

l = 0;

hold on;
for i = cellfun(@(x) s(x), {'O(0)', 'N(V)', 'N(-III)'})
    u = sol(:, :, i);
    
    if l
        u(u < 0) = 1e-7;
    end
    
    surf(u)
end
xlabel('depth');
ylabel('time');
zlabel('concentration');

if l
    set(gca, 'ZScale', 'log');
end

hold off;