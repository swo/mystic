s = species_map();
sol = run();

clf;

twod = 1;

hold on;
for i = cellfun(@(x) s(x), {'O', 'C', 'N+'})
    u = sol(:, :, i);
    
    surf(u)
end
xlabel('depth');
ylabel('time');
zlabel('concentration');

hold off;