s = species_map();
sol = run();

clf;

twod = 1;

hold on;
for i = cellfun(@(x) s(x), {'O(0)', 'N(-III)', 'N(V)', 'C(IV)', 'Fe(III)'})
    u = sol(:, :, i);
    
    surf(u)
end
xlabel('depth');
ylabel('time');
zlabel('concentration');

hold off;

if twod
    clf;
    x = linspace(0, 10, 5);
    f = @(x) sol(end,:, s(x));
    
    plot(x, f('O(0)'), x, f('N(-III)'), x, f('C(0)'), x, f('C(IV)'), x, f('photons'))
    legend('oxygen', 'nitrate', 'glucose', 'co2', 'photons')
end