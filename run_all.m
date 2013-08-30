s = species_map();
sol = run_mass();

clf;

twod = 1;

hold on;
for i = cellfun(@(x) s(x), {'O(0)', 'N(V)', 'N(-III)'})
    u = sol(:, :, i);
    
    surf(u)
end
xlabel('depth');
ylabel('time');
zlabel('concentration');

hold off;